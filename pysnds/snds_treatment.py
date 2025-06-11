'''
Class for detecting the presence of an event in the SNDS for a targeted population and determining its occurrence dates, including the first appearance.
'''

from .snds_query import SNDS_Query
import pandas as pd
import numpy as np

class SNDS_Treatment(SNDS_Query) :
    """
    Class for detecting the presence of an event in the SNDS for a targeted population and determining its occurrence dates.
    """

    def __init__(self, conn, df_ID_PATIENT):
        super().__init__(conn)

        if set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}:
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        self.df_ID_PATIENT = df_ID_PATIENT


    def Had_Treatment(self, dict_code, print_option=True) :
        '''
        Method to determine whether an event occurs for a targeted population.

        Parameters
        ----------
        dict_code : dict
            Dictionary of codes referring to the event of interest. Each key represents a code type 
            (possible keys: {'CCAM', 'CIP13', 'UCD', 'ICD10'}) and maps to a list of corresponding codes.
        print_option : bool, optional
            If True, prints the number of unique patients identified as having the event occur in each table (PMSI and DCIR) 
            and for each type of code considered. Default is True.

        Returns
        -------
        df_Treatment : DataFrame
            DataFrame containing the event response ('Response') for each patient in the targeted population ('BEN_IDT_ANO').
        '''

        if type(dict_code) != dict :
            raise ValueError("dict_code must be a dictionnary with keys 'CCAM', 'CIP13', 'UCD' and/or 'ICD10'.")

        unique_identifier = pd.DataFrame(columns=['BEN_IDT_ANO'])
        df_Treatment = self.df_ID_PATIENT.copy()
        df_Treatment['Response'] = 0

        for key in dict_code :

            if key == 'CCAM' :
                ccam_dcir = self.loc_ccam_dcir(self.df_ID_PATIENT, dict_code['CCAM'], print_option)
                ccam_pmsi = self.loc_ccam_pmsi(self.df_ID_PATIENT, dict_code['CCAM'], print_option)
            
                unique_identifier_CCAM = pd.DataFrame(pd.concat([ccam_dcir['BEN_IDT_ANO'], ccam_pmsi['BEN_IDT_ANO']], axis=0).drop_duplicates().reset_index(drop=True))
                unique_identifier = pd.DataFrame(pd.concat([unique_identifier['BEN_IDT_ANO'], unique_identifier_CCAM['BEN_IDT_ANO']], axis=0).drop_duplicates().reset_index(drop=True))
                
            if key == 'ICD10' :
                icd10_pmsi = self.loc_icd10_pmsi(self.df_ID_PATIENT, dict_code['ICD10'], print_option)
                unique_identifier = pd.DataFrame(pd.concat([unique_identifier['BEN_IDT_ANO'], icd10_pmsi['BEN_IDT_ANO']], axis=0).drop_duplicates().reset_index(drop=True))

            if key == 'UCD' :
                ucd_pmsi = self.loc_ucd_pmsi(self.df_ID_PATIENT, dict_code['UCD'], print_option)
                unique_identifier = pd.DataFrame(pd.concat([unique_identifier['BEN_IDT_ANO'], ucd_pmsi['BEN_IDT_ANO']], axis=0).drop_duplicates().reset_index(drop=True))

            if key == 'CIP13' :
                cip13_dcir = self.loc_cip_dcir(self.df_ID_PATIENT, dict_code['CIP13'], print_option)
                unique_identifier = pd.DataFrame(pd.concat([unique_identifier['BEN_IDT_ANO'], cip13_dcir['BEN_IDT_ANO']], axis=0).drop_duplicates().reset_index(drop=True))

        df_Treatment.loc[df_Treatment['BEN_IDT_ANO'].isin(unique_identifier['BEN_IDT_ANO']), 'Response'] = 1

        return df_Treatment


    def treatment_dates(self, dict_code) :
        '''
        Method to retrieve the occurrence dates of a specific event for the targeted population.

        Parameters
        ----------
        dict_code : dict
            Dictionary of codes referring to the event of interest. Each key represents a code type 
            (possible keys: {'CCAM', 'CIP13', 'UCD', 'ICD10'}) and maps to a list of corresponding codes.

        Returns
        -------
        df_date : DataFrame
            DataFrame containing the event based on the type of code ('COD_ACT', 'COD_DIAG', 'COD_UCD', 'COD_CIP'), 
            along with its occurrence date ('DATE'), for each patient in the targeted population ('BEN_IDT_ANO').
        '''

        if type(dict_code) != dict :
            raise ValueError("dict_code must be a dictionnary with keys 'CCAM', 'CIP13', 'UCD' and/or 'ICD10'.")

        df_date = pd.DataFrame(columns=['BEN_IDT_ANO', 'COD_ACT', 'COD_DIAG', 'COD_UCD', 'COD_CIP', 'DATE'])
        
        for key in dict_code :

            if key == 'CCAM' :
                ccam_dcir = self.loc_ccam_dcir(self.df_ID_PATIENT, dict_code['CCAM'], print_option=False)
                ccam_pmsi = self.loc_ccam_pmsi(self.df_ID_PATIENT, dict_code['CCAM'], print_option=False)
                ccam_pmsi['DATE'] = pd.to_datetime(ccam_pmsi['EXE_SOI_AMD']) + pd.to_timedelta(ccam_pmsi['ENT_DAT_DEL'], unit='days')
                
                df_ccam = pd.DataFrame({'BEN_IDT_ANO' : np.concatenate((ccam_dcir.BEN_IDT_ANO, ccam_pmsi.BEN_IDT_ANO)),
                                    'COD_ACT' : np.concatenate((ccam_dcir.CAM_PRS_IDE, ccam_pmsi.CDC_ACT)),
                                    'COD_DIAG' : np.nan,
                                    'COD_UCD' : np.nan,
                                    'COD_CIP' : np.nan,
                                    'DATE' : np.concatenate((pd.to_datetime(ccam_dcir.EXE_SOI_DTD).dt.strftime('%Y-%m-%d'), pd.to_datetime(ccam_pmsi.DATE).dt.strftime('%Y-%m-%d')))})
                
                df_date = pd.concat([df_date, df_ccam], ignore_index=True)            

            if key == 'ICD10' :
                icd10_pmsi = self.loc_icd10_pmsi(self.df_ID_PATIENT, dict_code['ICD10'], print_option=False)

                df_icd10 = pd.DataFrame({'BEN_IDT_ANO' : icd10_pmsi.BEN_IDT_ANO,
                                    'COD_ACT' : np.nan,
                                    'COD_DIAG' : icd10_pmsi.DGN_PAL,
                                    'COD_UCD' : np.nan,
                                    'COD_CIP' : np.nan,
                                    'DATE' : icd10_pmsi.EXE_SOI_AMD})
                
                df_date = pd.concat([df_date, df_icd10], ignore_index=True)
            
            if key == 'UCD' :
                ucd_pmsi = self.loc_ucd_pmsi(self.df_ID_PATIENT, dict_code['UCD'], print_option=False)
                ucd_pmsi['DATE'] = pd.to_datetime(ucd_pmsi['EXE_SOI_AMD']) + pd.to_timedelta(ucd_pmsi['DELAI'], unit='days')

                df_ucd = pd.DataFrame({'BEN_IDT_ANO' : ucd_pmsi.BEN_IDT_ANO,
                                    'COD_ACT' : np.nan,
                                    'COD_DIAG' : np.nan,
                                    'COD_UCD' : ucd_pmsi.UCD_COD,
                                    'COD_CIP' : np.nan,
                                    'DATE' : ucd_pmsi.EXE_SOI_AMD})
                
                df_date = pd.concat([df_date, df_ucd], ignore_index=True)

            if key == 'CIP13' :
                cip_dcir = self.loc_cip_dcir(self.df_ID_PATIENT, dict_code['CIP13'], print_option=False)

                df_cip = pd.DataFrame({'BEN_IDT_ANO' : cip_dcir.BEN_IDT_ANO,
                                    'COD_ACT' : np.nan,
                                    'COD_DIAG' : np.nan,
                                    'COD_UCD' : np.nan,
                                    'COD_CIP' : cip_dcir.PHA_CIP_C13,
                                    'DATE' : pd.to_datetime(cip_dcir.EXE_SOI_DTD).dt.strftime('%Y-%m-%d')})
                
                df_date = pd.concat([df_date, df_cip], ignore_index=True)

        return df_date
    

    def first_date_treatment(self, dict_code):
        '''
        Method to retrieve the firt occurrence date of a specific event for the targeted population.

        Parameters
        ----------
        dict_code : dict
            Dictionary of codes referring to the event of interest. Each key represents a code type 
            (possible keys: {'CCAM', 'CIP13', 'UCD', 'ICD10'}) and maps to a list of corresponding codes.

        Returns
        -------
        df_date : DataFrame
            DataFrame containing the first occurrence date of the targeted event ('DATE'), for each patient in the targeted population ('BEN_IDT_ANO').
        '''
        
        if type(dict_code) != dict :
            raise ValueError("dict_code must be a dictionnary with keys 'CCAM', 'CIP13', 'UCD' and/or 'ICD10'.")

        df_date = self.treatment_dates(dict_code)

        return df_date.groupby('BEN_IDT_ANO')['DATE'].min().reset_index()
