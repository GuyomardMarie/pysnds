from .snds_query import SNDS_Query
import pandas as pd
import numpy as np
from datetime import datetime

class SNDS_Treatment(SNDS_Query) :
    """
    Class for detecting the presence of an event in the SNDS for a targeted population and determining its occurrence dates.
    """

    def __init__(self, conn):
        super().__init__(conn)

    def Had_Treatment(self, dict_code, df_ID_PATIENT=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], print_option=True, dev=False) :
        '''
        Method to determine whether an event occurs for a targeted population.

        Parameters
        ----------
        dict_code : dict
            Dictionary of codes referring to the event of interest. Each key represents a code type 
            (possible keys: {'CCAM', 'CIP13', 'UCD', 'ATC', 'ICD10'}) and maps to a list of corresponding codes.
        df_ID_PATIENT : DataFrame
            DataFrame containing the columns "BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM", which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
        years : list
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        print_option : bool, optional
            If True, prints the number of unique patients identified as having the event occur in each table (PMSI and DCIR) 
            and for each type of code considered. Default is True.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.

        Returns
        -------
        df_Treatment : DataFrame
            DataFrame containing the event response ('Response') for each patient in the targeted population ('BEN_IDT_ANO').
        '''

        if (df_ID_PATIENT is not None) and (set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}):
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        if (years is None) or (not isinstance(years, list)):
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")
        
        if type(dict_code) != dict :
            raise ValueError("dict_code must be a dictionnary with keys 'CCAM', 'CIP13', 'UCD', 'ATC' and/or 'ICD10'.")

        unique_identifier = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA'])

        for key in dict_code :

            if key == 'CCAM' :
                ccam_dcir = self.loc_ccam_dcir(list_CCAM=dict_code['CCAM'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=print_option)
                ccam_pmsi = self.loc_ccam_pmsi(list_CCAM=dict_code['CCAM'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=print_option, dev=dev)
            
                unique_identifier_CCAM = pd.concat([ccam_dcir[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA']], ccam_pmsi[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA']]]).drop_duplicates().reset_index(drop=True)
                unique_identifier = pd.DataFrame(pd.concat([unique_identifier, unique_identifier_CCAM], axis=0).drop_duplicates().reset_index(drop=True))
                
            if key == 'ICD10' :
                icd10_pmsi = self.loc_icd10_pmsi(list_ICD10=dict_code['ICD10'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=print_option, dev=dev)
                unique_identifier_ICD10 = icd10_pmsi[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA']].drop_duplicates().reset_index(drop=True)
                unique_identifier = pd.DataFrame(pd.concat([unique_identifier, unique_identifier_ICD10], axis=0).drop_duplicates().reset_index(drop=True))
                
            if key == 'UCD' :
                ucd_pmsi = self.loc_ucd_pmsi(list_UCD=dict_code['UCD'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=print_option, dev=dev)
                ucd_dcir = self.loc_ucd_dcir(df_ID_PATIENT=df_ID_PATIENT, years=years, list_UCD=dict_code['UCD'], print_option=print_option)
                
                unique_identifier_UCD = pd.concat([ucd_dcir[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA']], ucd_pmsi[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA']]]).drop_duplicates().reset_index(drop=True)
                unique_identifier = pd.DataFrame(pd.concat([unique_identifier, unique_identifier_UCD], axis=0).drop_duplicates().reset_index(drop=True))

            if key == 'CIP13' :
                cip13_dcir = self.loc_cip_dcir(list_CIP13=dict_code['CIP13'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=print_option)
                unique_identifier_CIP13 = cip13_dcir[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA']].drop_duplicates().reset_index(drop=True)
                unique_identifier = pd.DataFrame(pd.concat([unique_identifier, unique_identifier_CIP13], axis=0).drop_duplicates().reset_index(drop=True))

            if key == 'ATC' :
                atc_pmsi = self.loc_atc_pmsi(list_ATC=dict_code['ATC'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=print_option, dev=dev)
                atc_dcir = self.loc_atc_dcir(df_ID_PATIENT=df_ID_PATIENT, years=years, list_ATC=dict_code['ATC'], print_option=print_option)
                
                unique_identifier_ATC = pd.concat([atc_dcir[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA']], atc_pmsi[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA']]]).drop_duplicates().reset_index(drop=True)
                unique_identifier = pd.DataFrame(pd.concat([unique_identifier, unique_identifier_ATC], axis=0).drop_duplicates().reset_index(drop=True))

                
        if df_ID_PATIENT is None :
            df_Treatment = unique_identifier.copy()
            df_Treatment["Response"] = 1
        else :
            df_Treatment = df_ID_PATIENT[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA']].copy()
            unique_identifier["BEN_RNG_GEM"] = unique_identifier["BEN_RNG_GEM"].astype(int) 
            df_Treatment["BEN_RNG_GEM"] = df_Treatment["BEN_RNG_GEM"].astype(int) 
            unique_identifier["Response"] = 1
            df_Treatment = df_Treatment.merge(unique_identifier, on=['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA'], how="left")
            df_Treatment["Response"] = df_Treatment["Response"].fillna(0).astype(int)

        if print_option==True :
            print(str(df_Treatment.Response.value_counts()[1]) + ' unique patients identified.')

        return df_Treatment


    def treatment_dates(self, dict_code, df_ID_PATIENT=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], dev=False) :
        '''
        Method to retrieve the occurrence dates of a specific event for the targeted population.

        Parameters
        ----------
        dict_code : dict
            Dictionary of codes referring to the event of interest. Each key represents a code type 
            (possible keys: {'CCAM', 'CIP13', 'UCD', 'ATC', 'ICD10'}) and maps to a list of corresponding codes.
        df_ID_PATIENT : DataFrame
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
        years : list
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.

        Returns
        -------
        df_date : DataFrame
            DataFrame containing the event based on the type of code ('COD_ACT', 'COD_DIAG', 'COD_UCD', 'COD_CIP', 'COD_ATC'), 
            along with its occurrence date ('DATE'), for each patient in the targeted population ('BEN_IDT_ANO').
        '''

        if (df_ID_PATIENT is not None) and (set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}):
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        if (years is None) or (not isinstance(years, list)):
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")
       
        if type(dict_code) != dict :
            raise ValueError("dict_code must be a dictionnary with keys 'CCAM', 'CIP13', 'UCD', 'ATC',  and/or 'ICD10'.")

        df_date = pd.DataFrame(columns=['BEN_IDT_ANO', 'COD_ACT', 'COD_DIAG', 'COD_UCD', 'COD_CIP', 'COD_ATC', 'DATE'])
        
        for key in dict_code :

            if key == 'CCAM' :
                ccam_dcir = self.loc_ccam_dcir(list_CCAM=dict_code['CCAM'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=False)
                ccam_pmsi = self.loc_ccam_pmsi(list_CCAM=dict_code['CCAM'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=False, dev=dev)
                ccam_pmsi['DATE'] = pd.to_datetime(ccam_pmsi['EXE_SOI_DTD']) #+ pd.to_timedelta(ccam_pmsi['ENT_DAT_DEL'], unit='days')
                
                df_ccam = pd.DataFrame({'BEN_IDT_ANO' : np.concatenate((ccam_dcir.BEN_IDT_ANO, ccam_pmsi.BEN_IDT_ANO)),
                                        'BEN_NIR_PSA' : np.concatenate((ccam_dcir.BEN_NIR_PSA, ccam_pmsi.BEN_NIR_PSA)),
                                        'BEN_RNG_GEM' : np.concatenate((ccam_dcir.BEN_RNG_GEM, ccam_pmsi.BEN_RNG_GEM)),
                                    'COD_ACT' : np.concatenate((ccam_dcir.CAM_PRS_IDE, ccam_pmsi.CDC_ACT)),
                                    'COD_DIAG' : np.nan,
                                    'COD_UCD' : np.nan,
                                    'COD_CIP' : np.nan,
                                    'COD_ATC' : np.nan,
                                    'DATE' : np.concatenate((pd.to_datetime(ccam_dcir.EXE_SOI_DTD).dt.strftime('%Y-%m-%d'), pd.to_datetime(ccam_pmsi.DATE).dt.strftime('%Y-%m-%d')))})
                
                df_date = pd.concat([d for d in [df_date, df_ccam] if not d.empty], ignore_index=True)            

            if key == 'ICD10' :
                icd10_pmsi = self.loc_icd10_pmsi(list_ICD10=dict_code['ICD10'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=False, dev=dev)

                df_icd10 = pd.DataFrame({'BEN_IDT_ANO' : icd10_pmsi.BEN_IDT_ANO,
                                         'BEN_NIR_PSA' : icd10_pmsi.BEN_NIR_PSA,
                                         'BEN_RNG_GEM' : icd10_pmsi.BEN_RNG_GEM,
                                    'COD_ACT' : np.nan,
                                    'COD_DIAG' : icd10_pmsi.DGN_PAL,
                                    'COD_UCD' : np.nan,
                                    'COD_CIP' : np.nan,
                                    'COD_ATC' : np.nan,
                                    'DATE' : pd.to_datetime(icd10_pmsi.EXE_SOI_DTD).dt.strftime('%Y-%m-%d')})
                
                df_date = pd.concat([d for d in [df_date, df_icd10] if not d.empty], ignore_index=True)
            
            if key == 'UCD' :
                ucd_pmsi = self.loc_ucd_pmsi(list_UCD=dict_code['UCD'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=False, dev=dev)
                ucd_pmsi['DATE'] = pd.to_datetime(ucd_pmsi['EXE_SOI_DTD']) #+ pd.to_timedelta(ucd_pmsi['DELAI'], unit='days')
                ucd_dcir = self.loc_ucd_dcir(list_UCD=dict_code['UCD'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=False)

                df_ucd = pd.DataFrame({'BEN_IDT_ANO' : np.concatenate((ucd_pmsi.BEN_IDT_ANO, ucd_dcir.BEN_IDT_ANO)),
                                       'BEN_NIR_PSA' : np.concatenate((ucd_pmsi.BEN_NIR_PSA, ucd_dcir.BEN_NIR_PSA)),
                                       'BEN_RNG_GEM' : np.concatenate((ucd_pmsi.BEN_RNG_GEM, ucd_dcir.BEN_RNG_GEM)),
                                    'COD_ACT' : np.nan,
                                    'COD_DIAG' : np.nan,
                                    'COD_UCD' : np.concatenate((ucd_pmsi.COD_UCD, ucd_dcir.COD_UCD)),
                                    'COD_CIP' : np.nan,
                                    'COD_ATC' : np.nan,
                                    'DATE' : np.concatenate((pd.to_datetime(ucd_pmsi.EXE_SOI_DTD).dt.strftime('%Y-%m-%d'), pd.to_datetime(ucd_dcir.EXE_SOI_DTD).dt.strftime('%Y-%m-%d')))})
                
                df_date = pd.concat([d for d in [df_date, df_ucd] if not d.empty], ignore_index=True)

            if key == 'CIP13' :
                cip_dcir = self.loc_cip_dcir(list_CIP13=dict_code['CIP13'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=False)

                df_cip = pd.DataFrame({'BEN_IDT_ANO' : cip_dcir.BEN_IDT_ANO,
                                       'BEN_NIR_PSA' : cip_dcir.BEN_NIR_PSA,
                                       'BEN_RNG_GEM' : cip_dcir.BEN_RNG_GEM,
                                    'COD_ACT' : np.nan,
                                    'COD_DIAG' : np.nan,
                                    'COD_UCD' : np.nan,
                                    'COD_CIP' : cip_dcir.PHA_CIP_C13,
                                    'COD_ATC' : np.nan,
                                    'DATE' : pd.to_datetime(cip_dcir.EXE_SOI_DTD).dt.strftime('%Y-%m-%d')})
                
                df_date = pd.concat([d for d in [df_date, df_cip] if not d.empty], ignore_index=True)
            
            if key == 'ATC' :
                atc_pmsi = self.loc_atc_pmsi(list_ATC=dict_code['ATC'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=False, dev=dev)
                atc_pmsi['DATE'] = pd.to_datetime(atc_pmsi['EXE_SOI_DTD']) #+ pd.to_timedelta(atc_pmsi['DELAI'], unit='days')
                atc_dcir = self.loc_atc_dcir(list_ATC=dict_code['ATC'], df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=False)

                if (atc_pmsi.shape[0]!=0) | (atc_dcir.shape[0]!=0) :
                    df_atc = pd.DataFrame({'BEN_IDT_ANO' : np.concatenate((atc_pmsi.BEN_IDT_ANO, atc_dcir.BEN_IDT_ANO)),
                                           'BEN_NIR_PSA' : np.concatenate((atc_pmsi.BEN_NIR_PSA, atc_dcir.BEN_NIR_PSA)),
                                           'BEN_RNG_GEM' : np.concatenate((atc_pmsi.BEN_RNG_GEM, atc_dcir.BEN_RNG_GEM)),
                                        'COD_ACT' : np.nan,
                                        'COD_DIAG' : np.nan,
                                        'COD_UCD' : np.nan,
                                        'COD_CIP' : np.nan,
                                        'COD_ATC' : np.concatenate((atc_pmsi.PHA_ATC_C07, atc_dcir.PHA_ATC_C07)),
                                        'DATE' : np.concatenate((pd.to_datetime(atc_pmsi.EXE_SOI_DTD).dt.strftime('%Y-%m-%d'), pd.to_datetime(atc_dcir.EXE_SOI_DTD).dt.strftime('%Y-%m-%d')))})
                    
                    df_date = pd.concat([d for d in [df_date, df_atc] if not d.empty], ignore_index=True)
                    df_date["DATE"] = pd.to_datetime(df_date["DATE"], format="%Y-%m-%d", errors="coerce")
                
        return df_date
    

    def first_date_treatment(self, dict_code, df_ID_PATIENT=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], dev=False):
        '''
        Method to retrieve the firt occurrence date of a specific event for the targeted population.

        Parameters
        ----------
        dict_code : dict
            Dictionary of codes referring to the event of interest. Each key represents a code type 
            (possible keys: {'CCAM', 'CIP13', 'UCD', 'ATC', 'ICD10'}) and maps to a list of corresponding codes.
        df_ID_PATIENT : DataFrame
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
        years : list
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.

        Returns
        -------
        df_date : DataFrame
            DataFrame containing the first occurrence date of the targeted event ('DATE'), for each patient in the targeted population ('BEN_IDT_ANO').
        '''
        
        if type(dict_code) != dict :
            raise ValueError("dict_code must be a dictionnary with keys 'CCAM', 'CIP13', 'UCD', 'ATC' and/or 'ICD10'.")

        if (df_ID_PATIENT is not None) and (set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}):
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        if (years is None) or (not isinstance(years, list)):
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")
       
        df_date = self.treatment_dates(dict_code=dict_code, df_ID_PATIENT=df_ID_PATIENT, years=years, dev=dev)
        df_date["DATE"] = pd.to_datetime(df_date["DATE"], format="%Y-%m-%d", errors="coerce")
        
        return df_date.groupby(['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'])['DATE'].min().reset_index()
