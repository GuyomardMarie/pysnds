from .snds_query import SNDS_Query
from .snds_treatment import SNDS_Treatment
import pkg_resources
import sqlite3
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from dateutil.relativedelta import relativedelta


class SNDS_BC(SNDS_Treatment) :
    """
    Class for computing Breast Cancer patients Statistical Analysis using the SNDS.
    """

    def __init__(self, conn, df_ID_PATIENT):

        if set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}:
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        self.conn = conn
        
        spark_session_type = None
        if "pyspark" in str(type(conn)): 
            try:
                from pyspark.sql import SparkSession
                spark_session_type = SparkSession
            except ImportError:
                raise ImportError("PySpark n'est pas installé, impossible d'utiliser Spark.")

        # Spark
        if spark_session_type is not None and isinstance(conn, spark_session_type):
            self.conn = conn
            self.backend = "spark"

        # Sqlite
        elif isinstance(conn, sqlite3.Connection):
            self.conn = conn
            self.backend = "sqlite"

        else:
            raise TypeError(
                f"Paramètre conn invalide : {type(conn)}. "
                f"Attendu SparkSession ou sqlite3.Connection."
            )
            
            
        self.df_ID_PATIENT = df_ID_PATIENT
        json_path = pkg_resources.resource_filename(__name__, 'BC_medical_codes.json')
        with open(json_path, 'r') as file:
            self.BC_medical_codes = json.load(file)

        self.SNDS_query = SNDS_Query(self.conn)
        self.SNDS_Treatment = SNDS_Treatment(self.conn)


    def Get_ID(self):
        '''
        Method for collecting unique identifiers of the population in IR_BEN_R.

        Returns
        -------
        df : DataFrame
            DataFrame containing the unique identifier of the population in a column named 'BEN_IDT_ANO'.
        '''

        query_ID_patient =  """
            SELECT DISTINCT
            A.BEN_IDT_ANO,
            A.BEN_NIR_PSA, 
            A.BEN_RNG_GEM
            FROM IR_BEN_R A

        """

        unique_id = self.GetQuery(query_ID_patient)
    
        print('We have ' + str(unique_id.shape[0]) + ' distinct identifiers, ie. patients in the database.')

        return unique_id
    
    
    
    def Date_Diag(self, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], dev=False):   
        '''
        Method for identifying the date of diagnosis.

        Parameters
        ----------
        years : list
        List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        dev : bool, optional
        If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.

        Returns
        -------
        df_dates_diag : DataFrame
        DataFrame containing the date of diagnosis (column Date_Diag) for each patient (BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_REM).
        '''
    
        if (years is None) or (not isinstance(years, list)):
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")

        # Init
        df_dates_diag = self.df_ID_PATIENT.copy()
        df_dates_diag["BEN_RNG_GEM"] = df_dates_diag["BEN_RNG_GEM"].astype(int)
        df_dates_diag['Date_Diag'] = np.nan

        # First Treatment date
        df_first_treatment = self.first_date_treatment(self.BC_medical_codes['All_BC_Codes'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, dev=dev)

        # Biopsy
        df_dates_diag_subset = df_dates_diag[df_dates_diag["Date_Diag"].isna()][['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM']].reset_index(drop=True)
        df_ccam_dcir_biopsy = self.loc_ccam_dcir(list_CCAM=self.BC_medical_codes['Diag_Proc']['Breast_core_biopsy']['CCAM'], df_ID_PATIENT=df_dates_diag_subset, years=years, print_option=False)
        df_ccam_pmsi_biopsy = self.loc_ccam_pmsi(list_CCAM=self.BC_medical_codes['Diag_Proc']['Breast_core_biopsy']['CCAM'], df_ID_PATIENT=df_dates_diag_subset, years=years, print_option=False, dev=dev)
        df_ccam_biopsy = pd.concat([df_ccam_dcir_biopsy[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'EXE_SOI_DTD']], df_ccam_pmsi_biopsy[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'EXE_SOI_DTD']]]).drop_duplicates().reset_index(drop=True)

        df_merged_biopsy = df_ccam_biopsy.merge(df_first_treatment, on=["BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"], suffixes=("_df_ccam_biopsy", "_df_first_treatment"))
        df_filtered_biopsy = df_merged_biopsy[
        (df_merged_biopsy["EXE_SOI_DTD"] > df_merged_biopsy["DATE"] - pd.DateOffset(years=1)) &
        (df_merged_biopsy["EXE_SOI_DTD"] < df_merged_biopsy["DATE"])
        ]
        df_result_biopsy = df_filtered_biopsy.groupby(["BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"], as_index=False).agg(DATE=("EXE_SOI_DTD","min"))
        df_result_biopsy["BEN_RNG_GEM"] = df_result_biopsy["BEN_RNG_GEM"].astype(int)

        df_merged = df_dates_diag.merge(df_result_biopsy, on=["BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"], how="left")
        df_merged["Date_Diag"] = df_merged["DATE"].combine_first(df_merged["Date_Diag"])
        df_dates_diag = df_merged.drop(columns="DATE")
        
        # Fine-needle aspiration cytology
        df_dates_diag_subset = df_dates_diag[df_dates_diag["Date_Diag"].isna()][['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM']].reset_index(drop=True)

        if df_dates_diag_subset.shape[0]!=0 :
            df_ccam_dcir_cytology = self.loc_ccam_dcir(list_CCAM=self.BC_medical_codes['Diag_Proc']['Fine_needle_aspiration_cytology']['CCAM'], df_ID_PATIENT=df_dates_diag_subset, years=years, print_option=False)
            df_ccam_pmsi_cytology = self.loc_ccam_pmsi(list_CCAM=self.BC_medical_codes['Diag_Proc']['Fine_needle_aspiration_cytology']['CCAM'], df_ID_PATIENT=df_dates_diag_subset, years=years, print_option=False, dev=dev)
            df_ccam_cytology = pd.concat([df_ccam_dcir_cytology[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'EXE_SOI_DTD']], df_ccam_pmsi_cytology[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'EXE_SOI_DTD']]]).drop_duplicates().reset_index(drop=True)

            if df_ccam_cytology.shape[0]!=0 :
                df_merged_cytology = df_ccam_cytology.merge(df_first_treatment, on=["BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"], suffixes=("_df_ccam_biopsy", "_df_first_treatment"))
                df_filtered_cytology = df_merged_cytology[
                (df_merged_cytology["EXE_SOI_DTD"] > df_merged_cytology["DATE"] - pd.DateOffset(years=1)) &
                (df_merged_cytology["EXE_SOI_DTD"] < df_merged_cytology["DATE"])
                ]
                df_result_cytology = df_filtered_cytology.groupby(["BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"], as_index=False).agg(DATE=("EXE_SOI_DTD","min"))
                df_result_cytology["BEN_RNG_GEM"] = df_result_cytology["BEN_RNG_GEM"].astype(int)
                
                df_merged = df_dates_diag.merge(df_result_cytology, on=["BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"], how="left")
                df_merged["Date_Diag"] = df_merged["DATE"].combine_first(df_merged["Date_Diag"])
                df_dates_diag = df_merged.drop(columns="DATE")
            
        # Breast Imaging Procedure
        df_dates_diag_subset = df_dates_diag[df_dates_diag["Date_Diag"].isna()][['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM']].reset_index(drop=True)

        if df_dates_diag_subset.shape[0]!=0 : 

            df_ccam_dcir_imagery = self.loc_ccam_dcir(list_CCAM=self.BC_medical_codes['Diag_Proc']['Breast_Imaging_Procedures']['All']['CCAM'], df_ID_PATIENT=df_dates_diag_subset, years=years, print_option=False)
            df_ccam_pmsi_imagery = self.loc_ccam_pmsi(list_CCAM=self.BC_medical_codes['Diag_Proc']['Breast_Imaging_Procedures']['All']['CCAM'], df_ID_PATIENT=df_dates_diag_subset, years=years, print_option=False, dev=dev)
            df_ccam_imagery = pd.concat([df_ccam_dcir_imagery[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'EXE_SOI_DTD']], df_ccam_pmsi_imagery[['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'EXE_SOI_DTD']]]).drop_duplicates().reset_index(drop=True)
            
            if df_ccam_imagery.shape[0]!=0 :
                df_ccam_imagery["EXE_SOI_DTD"] = pd.to_datetime(df_ccam_imagery["EXE_SOI_DTD"], errors="coerce")
                df_ccam_imagery["BEN_RNG_GEM"] = df_ccam_imagery["BEN_RNG_GEM"].astype(int)
                df_first_treatment["BEN_RNG_GEM"] = df_first_treatment["BEN_RNG_GEM"].astype(int)
                
                df_merged_imagery = df_ccam_imagery.merge(df_first_treatment, on=["BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"], suffixes=("_df_ccam_imagery", "_df_first_treatment"))
                df_filtered_imagery = df_merged_imagery[(df_merged_imagery["EXE_SOI_DTD"] < df_merged_imagery["DATE"])]

                def select_imagery_date_diag(group):
                    group = group.sort_values("EXE_SOI_DTD", ascending=False).reset_index(drop=True)
                    current = group.loc[0, "EXE_SOI_DTD"]

                    for d in group["EXE_SOI_DTD"][1:]:
                        if d > (current - relativedelta(months=1)):
                            current = d  
                        else:
                            break        

                    return current

                df_result_imagery = df_filtered_imagery.groupby(["BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"]).apply(select_imagery_date_diag).reset_index(name="DATE")
                df_result_imagery["BEN_RNG_GEM"] = df_result_imagery["BEN_RNG_GEM"].astype(int)
                
                df_merged = df_dates_diag.merge(df_result_imagery, on=["BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"], how="left")
                df_merged["Date_Diag"] = df_merged["DATE"].combine_first(df_merged["Date_Diag"])
                df_dates_diag = df_merged.drop(columns="DATE")

        # First treatment date
        df_dates_diag_subset = df_dates_diag[df_dates_diag["Date_Diag"].isna()][['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM']].reset_index(drop=True)

        if df_dates_diag_subset.shape[0]!=0 : 

            df_dates_diag_subset["BEN_RNG_GEM"] = df_dates_diag_subset["BEN_RNG_GEM"].astype(int)
            df_first_treatment["BEN_RNG_GEM"] = df_first_treatment["BEN_RNG_GEM"].astype(int)

            df_first_date = df_dates_diag_subset.merge(
            df_first_treatment[["BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM","DATE"]],
            on=["BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"],
            how="left"
            )
            df_first_date["BEN_RNG_GEM"] = df_first_date["BEN_RNG_GEM"].astype(int)
            
            df_merged = df_dates_diag.merge(df_first_date, on=["BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"], how="left")
            df_merged["Date_Diag"] = df_merged["DATE"].combine_first(df_merged["Date_Diag"])
            df_dates_diag = df_merged.drop(columns="DATE")
            
        return df_dates_diag


    def Get_AGE(self, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], dev=False):
        '''
        Method for collecting the population's age at the time of enrollment.

        Parameters
        ----------
        years : list
        List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        dev : bool, optional
        If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.

        Returns
        -------
        df_age : DataFrame
        DataFrame containing the age at the time of enrollment (column DATE_DIAG) for each patient (BEN_IDT_ANO).
        '''

        liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(self.df_ID_PATIENT.BEN_IDT_ANO))
        liste_bennirpsa_str = ', '.join(f"'{valeur}'" for valeur in np.unique(self.df_ID_PATIENT.BEN_NIR_PSA))

        if (years is None) or (not isinstance(years, list)):
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")

        if isinstance(years[0], datetime):
            flxmin_year = years[0].year
            year_deb = int(str(years[0].year)[-2:])
            deb = years[0]
        else:
            flxmin_year = years[0]
            year_deb = str(years[0])[-2:]
            deb = datetime(years[0], 1, 1)

        if isinstance(years[1], datetime):
            end = years[1]
            year_end = int(str(years[1].year)[-2:])
            flxmax_year = years[1].year
        else:
            end = datetime(years[1], 12, 31)
            year_end = str(years[1])[-2:]
            flxmax_year = years[1]

        flxmin = datetime(flxmin_year, 1, 1)
        flxmax_date = (end + relativedelta(months=6)).replace(day=1)
        vecflx = []
        current_date = flxmin
        while current_date <= flxmax_date:
            vecflx.append(current_date.strftime("%Y-%m-%d"))
            current_date += relativedelta(months=1)

        top_ER_PRS_F = True

        # if self.backend == 'sqlite':
        #     cursor = self.conn.cursor()
        #     cursor.execute(
        #     "SELECT name FROM sqlite_master WHERE type='table' AND name=?;",
        #     ('E_PRS_F',)
        #     )
        #     top_ER_PRS_F = cursor.fetchone() is not None

        if self.backend == 'spark':
            top_ER_PRS_F = self.conn.catalog.tableExists('ER_PRS_F')


        # DCIR
        age_dcir = pd.DataFrame(columns=['BEN_NIR_PSA', 'BEN_RNG_GEM', 'BEN_IDT_ANO', 'EXE_SOI_DTD', 'BEN_AMA_COD'])

        if top_ER_PRS_F==True :

            for flux in vecflx:
                query_AGE_DCIR = f"""

                    SELECT  
                    A.BEN_NIR_PSA,      
                    A.BEN_RNG_GEM,
                    B.BEN_IDT_ANO,  
                    A.EXE_SOI_DTD,
                    A.BEN_AMA_COD
                    FROM ER_PRS_F A

                    INNER JOIN IR_BEN_R B
                    ON B.BEN_NIR_PSA = A.BEN_NIR_PSA AND B.BEN_RNG_GEM = A.BEN_RNG_GEM

                    WHERE A.FLX_DIS_DTD = '{flux}' AND B.BEN_IDT_ANO IN ({liste_benidtano_str}) AND B.BEN_NIR_PSA IN ({liste_bennirpsa_str})
                """
                df_flux = self.GetQuery(query_AGE_DCIR)
                age_dcir = pd.concat([age_dcir, df_flux], ignore_index=True)
                age_dcir['EXE_SOI_DTD'] = pd.to_datetime(age_dcir['EXE_SOI_DTD']) 
                age_dcir.drop_duplicates(inplace=True)

        else :
            yy = list(range(int(flxmin_year), (int(flxmax_year)+1)))

            for year in yy :

                query_AGE_DCIR = f"""

                    SELECT 
                    A.BEN_NIR_PSA,      
                    A.BEN_RNG_GEM,
                    B.BEN_IDT_ANO,  
                    A.EXE_SOI_DTD,
                    A.BEN_AMA_COD
                    FROM ER_PRS_F_{year} A

                    INNER JOIN IR_BEN_R B
                    ON B.BEN_NIR_PSA = A.BEN_NIR_PSA AND B.BEN_RNG_GEM = A.BEN_RNG_GEM

                    WHERE B.BEN_IDT_ANO IN ({liste_benidtano_str}) AND B.BEN_NIR_PSA IN ({liste_bennirpsa_str})
                """
                df_flux = self.GetQuery(query_AGE_DCIR)
                age_dcir = pd.concat([age_dcir, df_flux], ignore_index=True)
                age_dcir['EXE_SOI_DTD'] = pd.to_datetime(age_dcir['EXE_SOI_DTD']) 
                age_dcir.drop_duplicates(inplace=True)

        # PMSI
        age_pmsi = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'EXE_SOI_DTD', 'AGE_ANN'])

        if dev==True:
            query_AGE_PMSI = f"""

                SELECT DISTINCT
                C.BEN_IDT_ANO, 
                C.BEN_RNG_GEM, 
                C.BEN_NIR_PSA, 
                A.EXE_SOI_DTD, 
                B.AGE_ANN
                FROM T_MCOaaC A

                INNER JOIN T_MCOaaB B
                ON A.ETA_NUM = B.ETA_NUM AND A.RSA_NUM = B.RSA_NUM

                INNER JOIN IR_BEN_R C
                ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                WHERE C.BEN_IDT_ANO IN ({liste_benidtano_str}) AND C.BEN_NIR_PSA IN ({liste_bennirpsa_str})
            """
            age_pmsi = self.GetQuery(query_AGE_PMSI)
            age_pmsi['EXE_SOI_DTD'] = pd.to_datetime(age_pmsi['EXE_SOI_DTD'])  
            age_pmsi.drop_duplicates(inplace=True)

        else :
            yy = list(range(int(year_deb), int(year_end)+1))

            for year in yy :
                table_nameB = f"T_MCO{year}B"
                table_nameC = f"T_MCO{year}C"
                query_AGE_PMSI = f"""

                    SELECT DISTINCT
                    C.BEN_IDT_ANO, 
                    C.BEN_RNG_GEM, 
                    C.BEN_NIR_PSA,  
                    A.EXE_SOI_DTD, 
                    B.AGE_ANN
                    FROM {table_nameC} A

                    INNER JOIN {table_nameB} B
                    ON A.ETA_NUM = B.ETA_NUM AND A.RSA_NUM = B.RSA_NUM

                    INNER JOIN IR_BEN_R C
                    ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                    WHERE C.BEN_IDT_ANO IN ({liste_benidtano_str}) AND C.BEN_NIR_PSA IN ({liste_bennirpsa_str})
                """
                df_flux = self.GetQuery(query_AGE_PMSI)
                age_pmsi = pd.concat([age_pmsi, df_flux], ignore_index=True)
                age_pmsi['EXE_SOI_DTD'] = pd.to_datetime(age_pmsi['EXE_SOI_DTD'])  
                age_pmsi.drop_duplicates(inplace=True)

        # AGE
        age_pmsi = age_pmsi.rename(columns={'AGE_ANN': 'AGE'})
        age_dcir = age_dcir.rename(columns={'BEN_AMA_COD': 'AGE'})

        cols = ['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM', 'EXE_SOI_DTD', 'AGE']
        df_age = pd.concat([age_dcir[cols], age_pmsi[cols]], ignore_index=True)
        df_age.drop_duplicates(inplace=True)
        df_age.sort_values(['BEN_IDT_ANO', 'AGE'])
        df_age.reset_index(drop=True, inplace=True)

        return df_age

    
    def Age_Diagnosis(self, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], dev=False):   
        '''
        Method for identifying the date of diagnosis and the age at this time.

        Parameters
        ----------
        years : list
        List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        dev : bool, optional
        If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.

        Returns
        -------
        df_diag : DataFrame
        DataFrame containing the date of diagnosis (column Date_Diag) and the age (AGE) for each patient (BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_REM).
        '''

        df_date_diag = self.Date_Diag(years=years, dev=dev)
        df_date_diag["BEN_RNG_GEM"] = df_date_diag["BEN_RNG_GEM"].astype(int)

        df_age_diag = self.Get_AGE(years=years, dev=dev)
        df_age_diag["BEN_RNG_GEM"] = df_age_diag["BEN_RNG_GEM"].astype(int)

        df_date_diag['Date_Diag'] = pd.to_datetime(df_date_diag['Date_Diag']).dt.date
        df_age_diag['EXE_SOI_DTD'] = pd.to_datetime(df_age_diag['EXE_SOI_DTD']).dt.date
        
        df_diag = df_age_diag.merge(df_date_diag, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')
        df_diag = df_diag[df_diag.EXE_SOI_DTD==df_diag.Date_Diag]
        df_diag.drop(['EXE_SOI_DTD'], axis=1, inplace=True)
        df_diag["AGE"] = df_diag["AGE"].astype(int)
        df_diag["Date_Diag"] = pd.to_datetime(df_diag["Date_Diag"], errors="coerce")
        
        return df_diag


    
    
    
    def treatment_setting(self, dict_treatment,  years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], dev=False):
        '''
        Method to determine the setting of a treatment : neoadjuvant (before the surgery) or adjuvant (after the surgery).

        Parameters
        ----------
        dict_treatment : dict
            Dictionary of codes referring to the treatment of interest. Each key represents a code type 
            (possible keys: {'CCAM', 'CIP13', 'UCD', 'ICD10'}) and maps to a list of corresponding codes.
        years : list, optional
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.

        Returns
        -------
        df_treatment_setting : DataFrame
            DataFrame containing the event response ('Setting') for each patient in the targeted population ('BEN_IDT_ANO'). 
            The values taken by the variable are 'No', 'Neoadjuvant' or 'Adjuvant'.
        '''
        
        if type(dict_treatment) != dict :
            raise ValueError("dict_treatment must be a dictionnary with keys 'CCAM', 'CIP13', 'UCD' and/or 'ICD10'.")

        # Compute first date of treatment
        surgery_date = self.first_date_treatment(dict_code=self.BC_medical_codes['Surgery_BC']['Surgery'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, dev=dev)
        treatment_date = self.first_date_treatment(dict_code=dict_treatment, df_ID_PATIENT=self.df_ID_PATIENT, years=years, dev=dev)

        merged = pd.merge(treatment_date, 
                        surgery_date, 
                        on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], 
                        how='left', 
                        suffixes=('_treatment', '_surgery'))

        # No Surgery
        merged['Setting'] = '0'
        merged.loc[merged['DATE_surgery'].isna(), 'Setting'] = 'Neoadjuvant'
        # Neoadjuvant
        merged.loc[np.where(merged['DATE_treatment'] <= merged['DATE_surgery'])[0], 'Setting'] = 'Neoadjuvant'
        # Adjuvant
        merged.loc[np.where(merged['DATE_treatment'] > merged['DATE_surgery'])[0], 'Setting'] = 'Adjuvant'

        df_treatment_setting = pd.merge(self.df_ID_PATIENT, merged, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left', suffixes=('', '_merged'))
        df_treatment_setting['Setting'] = df_treatment_setting['Setting'].fillna('No')
        
        return df_treatment_setting



    def Chemotherapy_Regimen(self, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], dev=False) :
        '''
        Method to determine the regimen of Chemotherapy

        Returns
        -------
        df_res : DataFrame
            DataFrame containing the event response ('Regimen') for each patient in the targeted population ('BEN_IDT_ANO'). 
            The values taken by the variable are 'No' (for No Chemotherapy), 'Unitherapy' and 'Bitherapy'. 
        years : list, optional
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.
        
        Returns
        -------
        df_res : DataFrame
            DataFrame containing the event response ('CT_Regimen') for each patient in the targeted population ('BEN_IDT_ANO'). 
            The values taken by the variable are 'No', 'Unitherapy' or 'Bitherapy'.
        '''

        df_chemotherapy = self.Had_Treatment(self.BC_medical_codes['CT'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, print_option=False, dev=dev)
        df_res = df_chemotherapy.copy()
        df_res['CT_Regimen'] = np.nan
        df_res["CT_Regimen"] = df_res["CT_Regimen"].astype(str)


        df_CT = self.treatment_dates(dict_code=self.BC_medical_codes['CT'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, dev=dev)
        df_CT.drop(['COD_CIP', 'COD_ATC'], axis=1, inplace=True)
        df_CT['DATE'] = pd.to_datetime(df_CT['DATE'])
        df_CT = df_CT.groupby(['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM', 'DATE'], as_index=False).agg({
            'COD_UCD': 'first', 
            'COD_ACT': 'first',
            'COD_DIAG': 'first'
        })
        df_CT = df_CT.drop_duplicates(subset=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM', 'DATE'])

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


        df_CT_Regimen = df_filtered.groupby(['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM']).apply(determine_CT_treatment).reset_index(name='Regimen')

        df_res = df_res.merge(df_CT_Regimen, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')
        df_res['CT_Regimen'].replace('nan', np.nan, inplace=True)
        df_res['CT_Regimen'] = df_res['CT_Regimen'].fillna(df_res['Regimen'])
        df_res = df_res.drop(columns=['Regimen'])

        return df_res



    def EndoctrineTherapy_Treatment(self, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], dev=False) : 
        '''
        Method to determine the regimen of Endoctrine Therapy.
        years : list, optional
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.
        
        Returns
        -------
        df_res : DataFrame
            DataFrame containing the event response ('ET_Regimen') for each patient in the targeted population ('BEN_IDT_ANO'). 
            The values taken by the variable are 'No' (for No Endoctrine Therapy), 'Tamoxifen', 'AI', 'Tamoxifen with Agonist', 'AI with Agonist', 'Tamoxifen followed by AI', 'AI followed by Tamoxifen'.
        '''

        df_endoctrine_therapy = self.Had_Treatment(self.BC_medical_codes['ET']['All'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, print_option=False, dev=dev)
        df_res = df_endoctrine_therapy.copy()
        df_res.loc[df_res['BEN_IDT_ANO'].isin(df_endoctrine_therapy[df_endoctrine_therapy.Response==0].BEN_IDT_ANO), 'ET_Regimen'] = 'No ET'

        df_ET = self.treatment_dates(self.BC_medical_codes['ET']['All'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, dev=dev)
        df_ET['DATE'] = pd.to_datetime(df_ET['DATE'])
        df_ET['COD_CIP'] = df_ET['COD_CIP'].astype(str)
        df_ET.replace('nan', np.nan, inplace=True)
        df_ET = df_ET[df_ET.COD_CIP.notna()] 
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


        #df_ET_Regimen = (df_ET.groupby('BEN_IDT_ANO').apply(determine_ET_treatment).reset_index().rename(columns={0: "Regimen"}))
        df_ET_Regimen = df_ET.groupby(['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM']).apply(determine_ET_treatment).reset_index(name='Regimen')
        df_res = df_res.merge(df_ET_Regimen, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='outer')
        df_res['Regimen'].replace('nan', np.nan, inplace=True)
        df_res['ET_Treatment'] = df_res['ET_Regimen'].fillna(df_res['Regimen'])
        df_res = df_res.drop(columns=['Regimen'])

        return df_res
    
    
    def BC_POP_Stat(self, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], dev=False) :
        '''
        Method to characterize the Breast Cancer Population.
        years : list, optional
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.
        
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
                - 'ET' : No : '0' / Yes : '1'
                - 'ET_Setting' : 'No', 'Neoadjuvant' or 'Adjuvant'
                - 'ET_Treatment' : 'No', 'Tamoxifen', 'AI', 'Tamoxifen with Agonist', 'AI with Agonist', 'Tamoxifen followed by AI', 'AI followed by Tamoxifen'
                - 'ET_Regimen' : 'No', 'Unitherapy' or 'Bitherapy'
        '''

        ### Age
        df = self.Age_Diagnosis(years=years, dev=dev)
        
        ### Nodal Status
        df_nodal_status = self.Had_Treatment(self.BC_medical_codes['Diag_NodalStatus'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, print_option=False, dev=dev)
        df_nodal_status.rename(columns={'Response' : 'Nodal_Status'}, inplace=True)
        df = df.merge(df_nodal_status, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')

        ### Surgery
        # Mastectomy
        df_masectomy = self.Had_Treatment(self.BC_medical_codes['Surgery_BC']['Mastectomy'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, print_option=False, dev=dev)
        df_masectomy.rename(columns={'Response' : 'Mastectomy'}, inplace=True)
        df = df.merge(df_masectomy, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')

        # Partial Mastectomy
        df_partial_masectomy = self.Had_Treatment(self.BC_medical_codes['Surgery_BC']['Partial_Mastectomy'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, print_option=False, dev=dev)
        df_partial_masectomy.rename(columns={'Response' : 'Partial_Mastectomy'}, inplace=True)
        df = df.merge(df_partial_masectomy, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')

        # Surgery
        df['Surgery'] = ((df['Mastectomy'] == 1) | (df['Partial_Mastectomy'] == 1)).astype(int)


        ### Chemotherapy
        # Yes / No
        df_CT = self.Had_Treatment(self.BC_medical_codes['CT'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, print_option=False, dev=dev)
        df_CT.rename(columns={'Response' : 'CT'}, inplace=True)
        df = df.merge(df_CT, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')

        # Setting
        df_CT_setting = self.treatment_setting(self.BC_medical_codes['CT'], years=years, dev=dev)[['BEN_IDT_ANO','BEN_NIR_PSA', 'BEN_RNG_GEM', 'Setting']]
        df_CT_setting.rename(columns={'Setting' : 'CT_Setting'}, inplace=True)
        df = df.merge(df_CT_setting, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')

        # Regimen
        df_CT_regimen = self.Chemotherapy_Regimen(years=years, dev=dev)[['BEN_IDT_ANO','BEN_NIR_PSA', 'BEN_RNG_GEM', 'CT_Regimen']]
        df = df.merge(df_CT_regimen, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')


        ### Radiotherapy
        # Yes / No
        df_RT = self.Had_Treatment(self.BC_medical_codes['RT'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, print_option=False, dev=dev)
        df_RT.rename(columns={'Response' : 'RT'}, inplace=True)
        df = df.merge(df_RT, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')

        # Setting
        df_RT_setting = self.treatment_setting(self.BC_medical_codes['RT'], years=years, dev=dev)[['BEN_IDT_ANO','BEN_NIR_PSA', 'BEN_RNG_GEM', 'Setting']]
        df_RT_setting.rename(columns={'Setting' : 'RT_Setting'}, inplace=True)
        df = df.merge(df_RT_setting, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')


        ### TT
        # Yes / No
        df_TT = self.Had_Treatment(self.BC_medical_codes['TT']['Pertuzumab'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, print_option=False, dev=dev)
        df_TT.rename(columns={'Response' : 'TT'}, inplace=True)
        df = df.merge(df_TT, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')

        # Setting
        df_TT_setting = self.treatment_setting(self.BC_medical_codes['TT']['Pertuzumab'], years=years, dev=dev)[['BEN_IDT_ANO','BEN_NIR_PSA', 'BEN_RNG_GEM', 'Setting']]
        df_TT_setting.rename(columns={'Setting' : 'TT_Setting'}, inplace=True)
        df = df.merge(df_TT_setting, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')

        
        ### ET
        # Yes / No
        df_ET = self.Had_Treatment(self.BC_medical_codes['ET']['All'], df_ID_PATIENT=self.df_ID_PATIENT, years=years, print_option=False, dev=dev)
        df_ET.rename(columns={'Response' : 'ET'}, inplace=True)
        df = df.merge(df_ET, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')

        # Setting
        df_ET_setting = self.treatment_setting(self.BC_medical_codes['ET']['All'], years=years, dev=dev)[['BEN_IDT_ANO','BEN_NIR_PSA', 'BEN_RNG_GEM', 'Setting']]
        df_ET_setting.rename(columns={'Setting' : 'ET_Setting'}, inplace=True)
        df = df.merge(df_ET_setting, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')

        # Treatment
        df_ET_treatment = self.EndoctrineTherapy_Treatment(years=years, dev=dev)
        df = df.merge(df_ET_treatment, on=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM'], how='left')

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
                return 0
        
        df_pathway = df_char.groupby(['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM']).apply(def_pathway).reset_index(name='Pathway')

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
            # Luminal
            if ((pathway.ET==1) & (pathway.TT==0)).all():
                return 'Luminal'
            # TNBC
            if ((pathway.CT==1) & (pathway.ET==0) & (pathway.TT==0)).all():
                return 'TNBC'
            # Unknown
            else :
                return 'Unknown'
        
        df_subtype = df_char.groupby(['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM']).apply(def_BC_type).reset_index(name='BC_SubType')

        return df_subtype


    def statistical_analyses(self, df_final=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], dev=False, save_option=True, pathway=True, age_range=True, path='') :
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
        years : list, optional
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.
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
            df_char = self.BC_POP_Stat(years=years, dev=dev)
            df_pathway = self.therapeutic_pathway(df_char)
            df_subtype = self.BC_subtype(df_char)
            df_final = df_char.set_index('ID_PATIENT').join([df_pathway.set_index('ID_PATIENT'), df_subtype.set_index('ID_PATIENT')], how='inner').reset_index()

        # General Statistics
        if (pathway==False) & (age_range==False) :

            #print('General Statistical Analysis')
            #print('----------------------------')

            general_stat = {}

            # Pathway
            general_stat['Pathway'] = round(df_final.Pathway.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)

            # BC Subtype
            general_stat['BC_SubType'] = round(df_final.BC_SubType.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)
            
            # Nodal Status
            general_stat['Nodal_Status'] = round(df_final.Nodal_Status.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)

            # Surgery
            general_stat['Mastectomy'] = round(df_final.Mastectomy.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)
            general_stat['Partial_Mastectomy'] = round(df_final.Partial_Mastectomy.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)

            # CT
            general_stat['CT'] = round(df_final.CT.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)
            general_stat['CT_Setting'] = round(df_final.CT_Setting.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)
            general_stat['CT_Regimen'] = round(df_final.CT_Regimen.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)

            # RT
            general_stat['RT'] = round(df_final.RT.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)
            general_stat['RT_Setting'] = round(df_final.RT_Setting.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)

            # TT
            general_stat['TT'] = round(df_final.TT.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)
            general_stat['TT_Setting'] = round(df_final.TT_Setting.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)

            # ET
            general_stat['ET'] = round(df_final.ET.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)
            general_stat['ET_Setting'] = round(df_final.ET_Setting.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)
            general_stat['ET_Treatment'] = round(df_final.ET_Treatment.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)
            general_stat['ET_Regimen'] = round(df_final.ET_Regimen.value_counts(normalize=True).sort_index(key=lambda x:x.map(lambda v: (0, v) if isinstance(v, int) else (1, str(v)))) * 100,2)

            if save_option == True:
                with pd.ExcelWriter(path+'general_stat.xlsx', engine='openpyxl') as writer:
                    for sheet_name, df in general_stat.items():
                        df.to_excel(writer, sheet_name=sheet_name, index=True)

            return general_stat


        # statistic by age ranges
        if (age_range == True) & (pathway == False) :

            stat_age = {}
            df_final['Age_range'] = pd.cut(df_final['AGE'], 
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
        if (age_range == False) & (pathway == True) :

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
            df_final['Age_range'] = pd.cut(df_final['AGE'], 
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

            return stat_pathway_age



    def vizualisation_pop(self, var_x, var_y, df_final=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], dev=False, ax=None): 
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
        years : list, optional
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.
        ax : ax subplots, optional
            Axis where to plot the figure in a subplot.
        '''

        # Final characterization of the Breast Cancer population
        if df_final is None :
            df_char = self.BC_POP_Stat(years=years, dev=dev)
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
                if percentage > 0:  
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
                if percentage > 0:  
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


