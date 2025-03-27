from .snds_query import SNDS_Query
from .snds_treatment import SNDS_Treatment
import pkg_resources
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class SNDS_BC(SNDS_Treatment) :
    """
    Class for computing Breast Cancer patients Statistical Analysis using the SNDS.
    """

    def __init__(self, conn, df_ID_PATIENT):
        self.conn = conn
        self.df_ID_PATIENT = df_ID_PATIENT
        json_path = pkg_resources.resource_filename(__name__, 'BC_medical_codes.json')
        with open(json_path, 'r') as file:
            self.BC_medical_codes = json.load(file)

        self.SNDS_query = SNDS_Query(self.conn)
        self.SNDS_Treatment = SNDS_Treatment(self.conn, self.df_ID_PATIENT)

    def treatment_setting(self, dict_code_treatment):
        '''
        Method to determine the setting of a treatment : neoadjuvant (before the surgery) or adjuvant (after the surgery).

        Parameters
        ----------
        dict_code_treatment : dict
            Dictionary of codes referring to the treatment of interest. Each key represents a code type 
            (possible keys: {'CCAM', 'CIP13', 'UCD', 'ICD10'}) and maps to a list of corresponding codes.

        Returns
        -------
        df_treatment_setting : DataFrame
            DataFrame containing the event response ('Setting') for each patient in the targeted population ('BEN_IDT_ANO'). 
            The values taken by the variable are 'No' (for no surgery), 'Neoadjuvant' or 'Adjuvant'.
        '''
        
        # Compute first date of treatment
        surgery_date = self.first_date_treatment(self.BC_medical_codes['Surgery_BC']['Surgery'])
        treatment_date = self.first_date_treatment(dict_code_treatment)

        merged = pd.merge(treatment_date[['BEN_IDT_ANO', 'DATE']], 
                        surgery_date[['BEN_IDT_ANO', 'DATE']], 
                        on='BEN_IDT_ANO', 
                        how='left', 
                        suffixes=('_treatment', '_surgery'))

        # No Surgery
        merged['Setting'] = 0
        merged.loc[merged['DATE_surgery'].isna(), 'Setting'] = 'Neoadjuvant'
        # Neoadjuvant
        merged.loc[np.where(merged['DATE_treatment'] <= merged['DATE_surgery'])[0], 'Setting'] = 'Neoadjuvant'
        # Adjuvant
        merged.loc[np.where(merged['DATE_treatment'] > merged['DATE_surgery'])[0], 'Setting'] = 'Adjuvant'

        df_treatment_setting = pd.DataFrame({'BEN_IDT_ANO': self.df_ID_PATIENT.BEN_IDT_ANO})
        df_treatment_setting = pd.merge(df_treatment_setting, merged[['BEN_IDT_ANO', 'Setting']], on='BEN_IDT_ANO', how='left', suffixes=('', '_merged'))
        df_treatment_setting['Setting'] = df_treatment_setting['Setting'].fillna('No')
        
        return df_treatment_setting



    def Chemotherapy_Regimen(self) :
        '''
        Method to determine the regimen of Chemotherapy

        Returns
        -------
        df_res : DataFrame
            DataFrame containing the event response ('Regimen') for each patient in the targeted population ('BEN_IDT_ANO'). 
            The values taken by the variable are 'No' (for No Chemotherapy), 'Unitherapy' and 'Bitherapy'. 
        '''

        df_chemotherapy = self.Had_Treatment(self.BC_medical_codes['CT'], print_option=False)
        df_res = df_chemotherapy.copy()
        df_CT = self.treatment_dates(self.BC_medical_codes['CT'])
        df_res['CT_Regimen'] = np.nan
        df_CT['DATE'] = pd.to_datetime(df_CT['DATE'])

        df_CT = df_CT.groupby(['BEN_IDT_ANO', 'DATE'], as_index=False).agg({
            'COD_UCD': 'first', 
            'COD_ACT': 'first',
            'COD_DIAG': 'first',
            'COD_CIP': 'first'
        })
        df_CT = df_CT.drop_duplicates(subset=['BEN_IDT_ANO', 'DATE'])

        # No CT
        df_res.loc[df_res['BEN_IDT_ANO'].isin(df_chemotherapy[df_chemotherapy.Response==0].BEN_IDT_ANO), 'CT_Regimen'] = 'No'

        # Only one CT
        df_res.loc[df_res['BEN_IDT_ANO'].isin(np.unique(df_CT.iloc[np.where(df_CT.groupby('BEN_IDT_ANO').size()==1)[0]].BEN_IDT_ANO)), 'CT_Regimen'] = 'One Treatment'

        # Patient with CT
        df_CT = df_CT.sort_values(by=['BEN_IDT_ANO', 'DATE'])
        df_CT['DATE'] = pd.to_datetime(df_CT['DATE'])
        df_CT['DIFF'] = df_CT.groupby('BEN_IDT_ANO')['DATE'].diff().dt.days
        df_CT['DIFF'] = df_CT['DIFF'].fillna(0).astype(int)

        # Round numbers : [6, 8] : 7, [13, 15] : 14 and [20, 22] : 21
        def round_diff(value):
            if value in [13, 15]:  
                return 14
            elif value in [20, 22]:  
                return 21
            elif value in [6, 8]: 
                return 7
            else:
                return value  

        df_CT['DIFF'] = df_CT['DIFF'].apply(round_diff)


        df_CT = df_CT.sort_values(by=['BEN_IDT_ANO', 'DATE'])
        df_filtered = df_CT[df_CT['DIFF'] != 0]


        def determine_CT_treatment(group):
            '''
            Method to determine the treatment of Chemotherapy.

            Parameters
            ----------
            group : Dataframe
                Identified Chemotherapy's codes for a patient, with corresponding dates.
            Returns
            -------
            'Unitherapy' or 'Bitherapy'. 
            '''

            diffs = group['DIFF'].tolist()

            if all(diff == 7 for diff in diffs):
                return 'Unitherapy'

            indices_of_14 = [i for i, diff in enumerate(diffs) if diff == 14]
            if len(indices_of_14) >= 2:
                remaining_diffs = [diff for i, diff in enumerate(diffs[indices_of_14[1] + 1:]) if diff != 14]
                if all(d in [7, 21] for d in remaining_diffs):
                    return 'Bitherapy'

            return 'Unknown'


        df_CT_Regimen = df_filtered.groupby('BEN_IDT_ANO').apply(determine_CT_treatment).reset_index(name='Regimen')

        df_res = df_res.merge(df_CT_Regimen[['BEN_IDT_ANO', 'Regimen']], on='BEN_IDT_ANO', how='left')
        df_res['CT_Regimen'] = df_res['CT_Regimen'].fillna(df_res['Regimen'])
        df_res = df_res.drop(columns=['Regimen'])

        return df_res



    def EndoctrineTherapy_Treatment(self) : 
        '''
        Method to determine the regimen of Endoctrine Therapy.

        Returns
        -------
        df_res : DataFrame
            DataFrame containing the event response ('Regimen') for each patient in the targeted population ('BEN_IDT_ANO'). 
            The values taken by the variable are 'No' (for No Endoctrine Therapy), 'Tamoxifen', 'AI', 'Tamoxifen with Agonist', 'AI with Agonist', 'Tamoxifen followed by AI', 'AI followed by Tamoxifen'.
        '''
    

        df_endoctrine_therapy = self.Had_Treatment(self.BC_medical_codes['ET']['All'], print_option=False)
        df_res = df_endoctrine_therapy.copy()
        df_res.loc[df_res['BEN_IDT_ANO'].isin(df_endoctrine_therapy[df_endoctrine_therapy.Response==0].BEN_IDT_ANO), 'ET_Regimen'] = 'No ET'

        df_ET = self.treatment_dates(self.BC_medical_codes['ET']['All'])
        df_ET['DATE'] = pd.to_datetime(df_ET['DATE'])
        df_ET['COD_CIP'] = df_ET['COD_CIP'].astype(str)
        df_ET = df_ET.sort_values(by=['BEN_IDT_ANO', 'DATE'])

        def determine_ET_treatment(group):
            '''
            Method to determine the treatment of Endoctrine Therapy.

            Parameters
            ----------
            group : Dataframe
                Identified Endoctrine Therapy's codes for a patient, with corresponding dates.
            Returns
            -------
            'Tamoxifen', 'AI', 'Tamoxifen with Agonist', 'AI with Agonist', 'Tamoxifen followed by AI', 'AI followed by Tamoxifen'.
            '''

            if all(group.COD_CIP.isin(self.BC_medical_codes['ET']['Tamoxifen']['CIP13'])) :
                return 'Tamoxifen'
            if all(group.COD_CIP.isin(self.BC_medical_codes['ET']['Aromatase_Inhibitor']['CIP13'])) :
                return 'AI'
            
            if all(group.COD_CIP.isin(self.BC_medical_codes['ET']['Tamoxifen']['CIP13']) | group.COD_CIP.isin(self.BC_medical_codes['ET']['GnRH_agonists']['CIP13'])) :
                return 'Tamoxifen with Agonist'
            
            if all(group.COD_CIP.isin(self.BC_medical_codes['ET']['Aromatase_Inhibitor']['CIP13']) | group.COD_CIP.isin(self.BC_medical_codes['ET']['GnRH_agonists']['CIP13'])) :
                return 'AI with Agonist'
            
            if any(group.COD_CIP.isin(self.BC_medical_codes['ET']['Tamoxifen']['CIP13'])) and any(group.COD_CIP.isin(self.BC_medical_codes['ET']['Aromatase_Inhibitor']['CIP13'])) :
                group = group.drop_duplicates(subset=['BEN_IDT_ANO', 'COD_CIP'], keep='first')
                if (group.COD_CIP.tolist()[0] in self.BC_medical_codes['ET']['Tamoxifen']['CIP13']) & (group.COD_CIP.tolist()[1] in self.BC_medical_codes['ET']['Aromatase_Inhibitor']['CIP13']) :
                    return 'Tamoxifen followed by AI'
                if (group.COD_CIP.tolist()[0] in self.BC_medical_codes['ET']['Aromatase_Inhibitor']['CIP13']) & (group.COD_CIP.tolist()[1] in self.BC_medical_codes['ET']['Tamoxifen']['CIP13']) :
                    return 'AI followed by Tamoxifen'
            
            else : 
                return 'Unknown'


        df_ET_Regimen = df_ET.groupby('BEN_IDT_ANO').apply(determine_ET_treatment).reset_index(name='Regimen')
        df_res = df_res.merge(df_ET_Regimen[['BEN_IDT_ANO', 'Regimen']], on='BEN_IDT_ANO', how='outer')
        df_res['ET_Treatment'] = df_res['ET_Regimen'].fillna(df_res['Regimen'])
        df_res = df_res.drop(columns=['Regimen'])

        return df_res
    
    
    def BC_POP_Stat(self) :
        '''
        Method to characterize the Breast Cancer Population.

        Returns
        -------
        df : DataFrame
            DataFrame containing cacacterization of the therapeutical pathway for each patient in the targeted population ('BEN_IDT_ANO'). 
            For each patient we have :
            - Age (at diagnosis, ie. at the first treatment) : 'Age'
            - Nodal Status 
                'Nodal_Status': '0' for Negative, '1' for Positive
            - Surgery
                - 'Mastectomy' : No : '0' / Yes : '1'
                - 'Partial_Mastectomy' : No : '0' / Yes : '1'
                - 'Surgery' : No : '0' / Yes : '1'
            - Chemotherapy  
                - 'CT' : No : '0' / Yes : '1'
                - 'CT_Setting' : 'No', 'Neoadjuvant' or 'Adjuvant'
                - 'CT_Regimen' : 'No', 'Unitherapy' or 'Bitherapy'. 
            - Radiotherapy
                - 'RT' : No : '0' / Yes : '1'
                - 'RT_Setting' : 'No', 'Neoadjuvant' or 'Adjuvant'
            - Targeted Therapy
                - 'TT' : No : '0' / Yes : '1'
                - 'TT_Setting' : 'No', 'Neoadjuvant' or 'Adjuvant'
            - Endoctrine Therapy
                - 'CT' : No : '0' / Yes : '1'
                - 'CT_Setting' : 'No', 'Neoadjuvant' or 'Adjuvant'
                - 'CT_Treatment' : 'No', 'Tamoxifen', 'AI', 'Tamoxifen with Agonist', 'AI with Agonist', 'Tamoxifen followed by AI', 'AI followed by Tamoxifen'
                 - 'CT_Regimen' : 'No', 'Unitherapy' or 'Bitherapy'
        '''

        df = pd.DataFrame({'ID_PATIENT' : self.df_ID_PATIENT['BEN_IDT_ANO']})

        ### Age
        df = pd.merge(df, self.Get_AGE(self.df_ID_PATIENT), left_on='ID_PATIENT', right_on='BEN_IDT_ANO', how='left')
        df.drop(columns=['BEN_IDT_ANO'], inplace=True)

        ### Nodal Status
        df['Nodal_Status'] = self.Had_Treatment(self.BC_medical_codes['Diag_NodalStatus'], print_option=False)['Response']


        ### Surgery
        # Mastectomy
        df['Mastectomy'] = self.Had_Treatment(self.BC_medical_codes['Surgery_BC']['Mastectomy'], print_option=False)['Response']

        # Partial Mastectomy
        df['Partial_Mastectomy'] = self.Had_Treatment(self.BC_medical_codes['Surgery_BC']['Partial_Mastectomy'], print_option=False)['Response']

        # Surgery
        df['Surgery'] = ((df['Mastectomy'] == 1) | (df['Partial_Mastectomy'] == 1)).astype(int)


        ### Chemotherapy
        # Yes / No
        df['CT'] = self.Had_Treatment(self.BC_medical_codes['CT'], print_option=False)['Response']

        # Setting
        df['CT_Setting'] = self.treatment_setting(self.BC_medical_codes['CT'])['Setting']

        # Regimen
        df['CT_Regimen'] = self.Chemotherapy_Regimen()['CT_Regimen']


        ### Radiotherapy
        # Yes / No
        df['RT'] = self.Had_Treatment(self.BC_medical_codes['RT'], print_option=False)['Response']

        # Setting
        df['RT_Setting'] = self.treatment_setting(self.BC_medical_codes['RT'])['Setting']


        ### TT
        # Yes / No
        df['TT'] = self.Had_Treatment(self.BC_medical_codes['TT']['Pertuzumab'], print_option=False)['Response']
        
        # Setting
        df['TT_Setting'] = self.treatment_setting(self.BC_medical_codes['TT']['Pertuzumab'])['Setting']

        
        ### ET
        # Yes / No
        df['ET'] = self.Had_Treatment(self.BC_medical_codes['ET']['All'], print_option=False)['Response']

        # Setting
        df['ET_Setting'] = self.treatment_setting(self.BC_medical_codes['ET']['All'])['Setting']

        # Treatment
        df['ET_Treatment'] = self.EndoctrineTherapy_Treatment()['ET_Treatment']

        # Regimen
        df['ET_Regimen'] = np.select([df['ET_Treatment'].isin(['AI', 'Tamoxifen']),
                                      df['ET_Treatment'].isin(['AI followed by Tamoxifen', 'AI with Agonist', 'Tamoxifen followed by AI', 'Tamoxifen with Agonist']),
                                      df['ET_Treatment'].isin(['No ET'])], ['Unitherapy', 'Bitherapy', 'No ET'], default=np.nan)

        return df




    def therapeutic_pathway(self, df_char):
        '''
        Method to determine the therapeutical pathway for each breast cancer patient.

        Parameters
        ----------
        df_char : Dataframe
            Characterization of the population resulting from BC_POP_Stat().

        Returns
        -------
        df_pathway : Dataframe
            Therapeutical Pathway for each patient : 1, 2, ..., 10, 'Unknown'.
        '''

        def def_pathway(pathway):
            # Without CT
            if ((pathway.CT==0) & (pathway.RT==0) & (pathway.TT==0) & (pathway.ET==0)).all():
                return 1
            if ((pathway.CT==0) & (pathway.RT==0) & (pathway.TT==0) & (pathway.ET==1)).all():
                return 2
            if ((pathway.CT==0) & (pathway.RT==1) & (pathway.TT==0) & (pathway.ET==0)).all():
                return 3
            if ((pathway.CT==0) & (pathway.RT==1) & (pathway.TT==0) & (pathway.ET==1)).all():
                return 4

            # Neoadjuvant CT
            if ((pathway.CT==1) & (pathway.CT_Setting=='Neoadjuvant') & (pathway.RT==1) & (pathway.TT==0) & (pathway.ET==1)).all():
                return 5
            if ((pathway.CT==1) & (pathway.CT_Setting=='Neoadjuvant') & (pathway.RT==1) & (pathway.TT==0) & (pathway.ET==0)).all():
                return 6
            if ((pathway.CT==1) & (pathway.CT_Setting=='Neoadjuvant') & (pathway.RT==1) & (pathway.TT==1) & (pathway.TT_Setting=='Neoadjuvant') & (pathway.ET==1)).all():
                return 7

            # Adjuvant CT
            if ((pathway.CT==1) & (pathway.CT_Setting=='Adjuvant') & (pathway.RT==1) & (pathway.TT==1) & (pathway.TT_Setting=='Adjuvant') & (pathway.ET==1)).all():
                return 8
            if ((pathway.CT==1) & (pathway.CT_Setting=='Adjuvant') & (pathway.RT==1) & (pathway.TT==0) & (pathway.ET==0)).all():
                return 9
            if ((pathway.CT==1) & (pathway.CT_Setting=='Adjuvant') & (pathway.RT==1) & (pathway.TT==0) & (pathway.ET==1)).all():
                return 10
            
            else :
                return 'Unknown'
        
        df_pathway = df_char.groupby('ID_PATIENT').apply(def_pathway).reset_index(name='Pathway')

        return df_pathway
    

    def BC_subtype(self, df_char):
        '''
        Method to determine the Breast Cancer Subtype for each patient.

        Parameters
        ----------
        df_char : Dataframe
            Characterization of the population resulting from BC_POP_Stat().

        Returns
        -------
        df_subtype : Dataframe
            Subtype of Breast Cancer for each patient : 'HER2', 'Luminal', 'TNBC' or' 'Unknown'.
        '''

        def def_BC_type(pathway):
            # HER2
            if ((pathway.TT==1)).all():
                return 'HER2'
            if ((pathway.ET==1) & (pathway.TT==0)).all():
                return 'Luminal'
            if ((pathway.CT==1) & (pathway.ET==0) & (pathway.TT==0)).all():
                return 'TNBC'
            else :
                return 'Unknown'
        
        df_subtype = df_char.groupby('ID_PATIENT').apply(def_BC_type).reset_index(name='BC_SubType')

        return df_subtype


    def statistical_analyses(self, df_final = None, save_option=True, pathway=True, age_range=True, path='') :
        '''
        Method to make a statistical analysis of the breast cancer population.

        Parameters
        ----------
        df_final : Dataframe, optional
            Final Characterization of the population. Concatenation of :
                - df_char resulting from BC_POP_Stat()
                - df_pathway resulting from therapeutic_pathway()
                - df_subtype resulting from BC_subtype()
            If None it will be done directly in this method. Default is None.
        save_option : bool, optional
            If True, saves the resulting dataframe in excel format. Default is True.
        pathway : bool, optional
            If True, makes statistics regarding the pathway distribution in the population. Default is True.
        age_range : bool, optional
            If True, makes statistics regarding the age ange distribution (<50, [50-60[, [60, 70[, >=70) in the population. Default is True.
        path : str, optional
            Pathway where to save the resulting excel doc. Default is the current pathway.
        Returns
        -------
        Dataframe
            Statistics of the breast cancer population, regarding age range and pathways.
            - general_stat : for pathway = False, age_range = False.
                General Statistics
            - stat_age : for pathway = False, age_range = True.
                Statistic of the population regarding the age range.
            - stat_pathway : for pathway = True, age_range = False.
                Statistic of the population regarding the therapeutical pathway.
            - stat_age_pathway : for pathway = True, age_range = True.
                Statistic of the population regarding the age range and the therapeutical pathway.
            If save_option = True, the resulting dataframe is also saved in excel format in the pathway 'path'.
        '''
        
        # Final characterization of the Breast Cancer population
        if df_final is None :
            df_char = self.BC_POP_Stat()
            df_pathway = self.therapeutic_pathway(df_char)
            df_subtype = self.BC_subtype(df_char)
            df_final = df_char.set_index('ID_PATIENT').join([df_pathway.set_index('ID_PATIENT'), df_subtype.set_index('ID_PATIENT')], how='inner').reset_index()

        # General Statistics
        if (pathway==False) & (age_range==False) :

            #print('General Statistical Analysis')
            #print('----------------------------')

            general_stat = {}

            # Pathway
            general_stat['Pathway'] = round(df_final.Pathway.value_counts(normalize=True).sort_index() * 100,2)

            # BC Subtype
            general_stat['BC_SubType'] = round(df_final.BC_SubType.value_counts(normalize=True).sort_index() * 100,2)
            
            # Nodal Status
            general_stat['Nodal_Status'] = round(df_final.Nodal_Status.value_counts(normalize=True).sort_index() * 100,2)

            # Surgery
            general_stat['Mastectomy'] = round(df_final.Mastectomy.value_counts(normalize=True).sort_index() * 100,2)
            general_stat['Partial_Mastectomy'] = round(df_final.Partial_Mastectomy.value_counts(normalize=True).sort_index() * 100,2)

            # CT
            general_stat['CT'] = round(df_final.CT.value_counts(normalize=True).sort_index() * 100,2)
            general_stat['CT_Setting'] = round(df_final.CT_Setting.value_counts(normalize=True).sort_index() * 100,2)
            general_stat['CT_Regimen'] = round(df_final.CT_Regimen.value_counts(normalize=True).sort_index() * 100,2)

            # RT
            general_stat['RT'] = round(df_final.RT.value_counts(normalize=True).sort_index() * 100,2)
            general_stat['RT_Setting'] = round(df_final.RT_Setting.value_counts(normalize=True).sort_index() * 100,2)

            # TT
            general_stat['TT'] = round(df_final.TT.value_counts(normalize=True).sort_index() * 100,2)
            general_stat['TT_Setting'] = round(df_final.TT_Setting.value_counts(normalize=True).sort_index() * 100,2)

            # ET
            general_stat['ET'] = round(df_final.ET.value_counts(normalize=True).sort_index() * 100,2)
            general_stat['ET_Setting'] = round(df_final.ET_Setting.value_counts(normalize=True).sort_index() * 100,2)
            general_stat['ET_Treatment'] = round(df_final.ET_Treatment.value_counts(normalize=True).sort_index() * 100,2)
            general_stat['ET_Regimen'] = round(df_final.ET_Regimen.value_counts(normalize=True).sort_index() * 100,2)

            if save_option == True:
                with pd.ExcelWriter(path+'general_stat.xlsx', engine='openpyxl') as writer:
                    for sheet_name, df in general_stat.items():
                        df.to_excel(writer, sheet_name=sheet_name, index=True)

            return general_stat


        # statistic by age ranges
        if (age_range == True) & (pathway == False) :

            stat_age = {}
            df_final['Age_range'] = pd.cut(df_final['AGE_DIAG'], 
                         bins=[0, 50, 60, 70, float('inf')], 
                         labels=["<50", "[50-60[", "[60-70[", ">=70"], 
                         right=False)
            
            # Nodal Status
            stat_age['Nodal_Status'] = round(pd.crosstab(df_final['Nodal_Status'], df_final['Age_range'], margins=True, normalize=True) * 100,2)

            # BC Subtype
            stat_age['BC_SubType'] = round(pd.crosstab(df_final['BC_SubType'], df_final['Age_range'], margins=True, normalize=True) * 100,2)

            # Surgery
            stat_age['Mastectomy'] = round(pd.crosstab(df_final['Mastectomy'], df_final['Age_range'], margins=True, normalize=True) * 100,2)
            stat_age['Partial_Mastectomy'] = round(pd.crosstab(df_final['Partial_Mastectomy'], df_final['Age_range'], margins=True, normalize=True) * 100,2)

            # CT
            stat_age['CT'] = round(pd.crosstab(df_final['CT'], df_final['Age_range'], margins=True, normalize=True) * 100,2)
            stat_age['CT_Setting'] = round(pd.crosstab(df_final['CT_Setting'], df_final['Age_range'], margins=True, normalize=True) * 100,2)
            stat_age['CT_Regimen'] = round(pd.crosstab(df_final['CT_Regimen'], df_final['Age_range'], margins=True, normalize=True) * 100,2)

            # RT
            stat_age['RT'] = round(pd.crosstab(df_final['RT'], df_final['Age_range'], margins=True, normalize=True) * 100,2)
            stat_age['RT_Setting'] = round(pd.crosstab(df_final['RT_Setting'], df_final['Age_range'], margins=True, normalize=True) * 100,2)

            # TT
            stat_age['TT'] = round(pd.crosstab(df_final['TT'], df_final['Age_range'], margins=True, normalize=True) * 100,2)
            stat_age['TT_Setting'] = round(pd.crosstab(df_final['TT_Setting'], df_final['Age_range'], margins=True, normalize=True) * 100,2)

            # ET
            stat_age['ET'] = round(pd.crosstab(df_final['ET'], df_final['Age_range'], margins=True, normalize=True) * 100,2)
            stat_age['ET_Setting'] = round(pd.crosstab(df_final['ET_Setting'], df_final['Age_range'], margins=True, normalize=True) * 100,2)
            stat_age['ET_Treatment'] = round(pd.crosstab(df_final['ET_Treatment'], df_final['Age_range'], margins=True, normalize=True) * 100,2)
            stat_age['ET_Regimen'] = round(pd.crosstab(df_final['ET_Regimen'], df_final['Age_range'], margins=True, normalize=True) * 100,2)

            if save_option == True :
                with pd.ExcelWriter(path+'stat_by_age.xlsx', engine='openpyxl') as writer:
                    for sheet_name, df in stat_age.items():
                        df.to_excel(writer, sheet_name=sheet_name, index=True)
            
            return stat_age
        

        # statistic by pathway
        if (age_range == True) & (pathway == False) :

            stat_pathway = {}

            # Nodal Status
            stat_pathway['Nodal_Status'] = round(pd.crosstab(df_final['Nodal_Status'], df_final['Pathway'], margins=True, normalize=True) * 100,2)

            # BC Subtype
            stat_pathway['BC_SubType'] = round(pd.crosstab(df_final['BC_SubType'], df_final['Pathway'], margins=True, normalize=True) * 100,2)

            # Surgery
            stat_pathway['Mastectomy'] = round(pd.crosstab(df_final['Mastectomy'], df_final['Pathway'], margins=True, normalize=True) * 100,2)
            stat_pathway['Partial_Mastectomy'] = round(pd.crosstab(df_final['Partial_Mastectomy'], df_final['Pathway'], margins=True, normalize=True) * 100,2)

            # CT
            stat_pathway['CT'] = round(pd.crosstab(df_final['CT'], df_final['Pathway'], margins=True, normalize=True) * 100,2)
            stat_pathway['CT_Setting'] = round(pd.crosstab(df_final['CT_Setting'], df_final['Pathway'], margins=True, normalize=True) * 100,2)
            stat_pathway['CT_Regimen'] = round(pd.crosstab(df_final['CT_Regimen'], df_final['Pathway'], margins=True, normalize=True) * 100,2)

            # RT
            stat_pathway['RT'] = round(pd.crosstab(df_final['RT'], df_final['Pathway'], margins=True, normalize=True) * 100,2)
            stat_pathway['RT_Setting'] = round(pd.crosstab(df_final['RT_Setting'], df_final['Pathway'], margins=True, normalize=True) * 100,2)

            # TT
            stat_pathway['TT'] = round(pd.crosstab(df_final['TT'], df_final['Pathway'], margins=True, normalize=True) * 100,2)
            stat_pathway['TT_Setting'] = round(pd.crosstab(df_final['TT_Setting'], df_final['Pathway'], margins=True, normalize=True) * 100,2)

            # ET
            stat_pathway['ET'] = round(pd.crosstab(df_final['ET'], df_final['Pathway'], margins=True, normalize=True) * 100,2)
            stat_pathway['ET_Setting'] = round(pd.crosstab(df_final['ET_Setting'], df_final['Pathway'], margins=True, normalize=True) * 100,2)
            stat_pathway['ET_Treatment'] = round(pd.crosstab(df_final['ET_Treatment'], df_final['Pathway'], margins=True, normalize=True) * 100,2)
            stat_pathway['ET_Regimen'] = round(pd.crosstab(df_final['ET_Regimen'], df_final['Pathway'], margins=True, normalize=True) * 100,2)

            if save_option == True :
                with pd.ExcelWriter(path+'stat_by_pathway.xlsx', engine='openpyxl') as writer:
                    for sheet_name, df in stat_pathway.items():
                        df.to_excel(writer, sheet_name=sheet_name, index=True)
            
            return stat_pathway
        

        # Statistics by pathway and age range
        if (age_range == True) & (pathway == True) :
            
            stat_pathway_age = {}
            df_final['Age_range'] = pd.cut(df_final['AGE_DIAG'], 
                         bins=[0, 50, 60, 70, float('inf')], 
                         labels=["<50", "[50-60[", "[60-70[", ">=70"], 
                         right=False)
            
            # Nodal Status
            df = round(pd.crosstab(index=[df_final['Pathway'], df_final['Age_range']], columns=df_final['Nodal_Status'], margins=True, normalize=True) * 100,2)
            stat_pathway_age['Nodal_Status'] = {str(pathway): df.xs(pathway, level='Pathway') for pathway in df.index.get_level_values('Pathway').unique()}

            # Mastectomy
            df = round(pd.crosstab(index=[df_final['Pathway'], df_final['Age_range']], columns=df_final['Mastectomy'], margins=True, normalize=True) * 100,2)
            stat_pathway_age['Mastectomy'] = {str(pathway): df.xs(pathway, level='Pathway') for pathway in df.index.get_level_values('Pathway').unique()}

            # Partial Mastectomy
            df = round(pd.crosstab(index=[df_final['Pathway'], df_final['Age_range']], columns=df_final['Partial_Mastectomy'], margins=True, normalize=True) * 100,2)
            stat_pathway_age['Partial_Mastectomy'] = {str(pathway): df.xs(pathway, level='Pathway') for pathway in df.index.get_level_values('Pathway').unique()}

            # CT Regimen
            df = round(pd.crosstab(index=[df_final['Pathway'], df_final['Age_range']], columns=df_final['CT_Regimen'], margins=True, normalize=True) * 100,2)
            stat_pathway_age['CT_Regimen'] = {str(pathway): df.xs(pathway, level='Pathway') for pathway in df.index.get_level_values('Pathway').unique()}

            # ET Treatment
            df = round(pd.crosstab(index=[df_final['Pathway'], df_final['Age_range']], columns=df_final['ET_Treatment'], margins=True, normalize=True) * 100,2)
            stat_pathway_age['ET_Treatment'] = {str(pathway): df.xs(pathway, level='Pathway') for pathway in df.index.get_level_values('Pathway').unique()}

            # ET Regimen
            df = round(pd.crosstab(index=[df_final['Pathway'], df_final['Age_range']], columns=df_final['ET_Regimen'], margins=True, normalize=True) * 100,2)
            stat_pathway_age['ET_Regimen'] = {str(pathway): df.xs(pathway, level='Pathway') for pathway in df.index.get_level_values('Pathway').unique()}


            with pd.ExcelWriter("stat_pathway_age.xlsx") as writer:
                for pathway, sub_dict in stat_pathway_age.items():
                    for category, sub_df in sub_dict.items():
                        sheet_name = f"Pathway_{pathway}_{category}"
                        sub_df.to_excel(writer, sheet_name=sheet_name)




    def vizualisation_pop(self, var_x, var_y, df_final=None, ax=None): 
        '''
        Method to determine the visualize the distributions of the breast cancer population.

        Parameters
        ----------
        df_final : Dataframe, optional
            Final Characterization of the population. Concatenation of :
                - df_char resulting from BC_POP_Stat()
                - df_pathway resulting from therapeutic_pathway()
                - df_subtype resulting from BC_subtype()
            If None it will be done directly in this method. Default is None.
        var_x : str
            Name of the column in df_final for x axis.
        var_y : str
            Name of the column in df_final for y axis.
        ax : ax subplots, optional
            Axis where to plot the figure in a subplot.
        '''

        # Final characterization of the Breast Cancer population
        if df_final is None :
            df_char = self.BC_POP_Stat()
            df_pathway = self.therapeutic_pathway(df_char)
            df_subtype = self.BC_subtype(df_char)
            df_final = df_char.set_index('ID_PATIENT').join([df_pathway.set_index('ID_PATIENT'), df_subtype.set_index('ID_PATIENT')], how='inner').reset_index()

        # Figure
        unique_values = df_final[str(var_y)].unique()
        colors_list = ['gold', 'limegreen', 'cornflowerblue', 'mediumpurple', 'mediumvioletred', 'orange', 'peru'] if len(unique_values) > 2 else ['gold', 'mediumpurple']
        color_dict = {val: color for val, color in zip(unique_values, colors_list)}

        standalone = ax is None
        if standalone:
            fig, ax = plt.subplots(figsize=(5, 5))

        if var_x is None:
            counts = df_final[str(var_y)].value_counts(normalize=True) * 100
            bottom = 0
            for label, count in counts.items():
                ax.bar(str(var_y), count, bottom=bottom, color=color_dict[label], label=label)
                bottom += count

            for p, label in zip(ax.patches, counts.index):
                percentage = p.get_height()
                if percentage > 0:  # Afficher uniquement si > 0
                    x = p.get_x() + p.get_width() / 2
                    y = p.get_y() + p.get_height() / 2
                    ax.text(x, y, f'{percentage:.1f}%', ha='center', va='center', fontsize=12, color='black')

            ax.set_xticks([0])
            ax.set_xticklabels(["All"], fontsize=12)
            ax.set_ylim(0, 100)

        else:
            tab = pd.crosstab(df_final[var_x], df_final[str(var_y)], normalize='index') * 100
            tab.plot(kind='bar', stacked=True, ax=ax, color=[color_dict[col] for col in tab.columns])

            for p in ax.patches:
                percentage = p.get_height()
                if percentage > 0:  # Afficher uniquement si > 0
                    x = p.get_x() + p.get_width() / 2
                    y = p.get_y() + p.get_height() / 2
                    ax.text(x, y, f'{percentage:.1f}%', ha='center', va='center', fontsize=12)

            ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=12)
            ax.set_xlabel("")
            ax.set_ylim(0, 100)

        ax.legend().remove()

        if standalone:
            handles, labels = ax.get_legend_handles_labels()
            plt.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.1), ncol=len(labels), fontsize=12)
            plt.tight_layout()
            plt.show()
