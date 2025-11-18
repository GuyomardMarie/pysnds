import sqlite3
import pandas as pd
import numpy as np
from datetime import datetime
from dateutil.relativedelta import relativedelta


class SNDS_Query() :
    """
    Class for navigating and identifying population in the SNDS.    
    """

    def __init__(self, conn):
        '''
        Parameters
        ----------
        conn : sqlite3.Connection or pyspark.sql.session.SparkSession.
            Connection to the database.
        '''

        super(SNDS_Query, self).__init__()
        
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

        # Vérifier sqlite
        elif isinstance(conn, sqlite3.Connection):
            self.conn = conn
            self.backend = "sqlite"

        else:
            raise TypeError(
                f"Paramètre conn invalide : {type(conn)}. "
                f"Attendu SparkSession ou sqlite3.Connection."
            )
            
            

    def GetQuery(self, query):
        '''
        Method to execute a SQL query.
        
        Parameters
        ----------
        query : string
            SQL query.

        Returns
        -------
        df : DataFrame
            DataFrame containing the targeted population.
        '''
        
        if self.backend == "spark":
            df = self.conn.sql(query)
            return df.toPandas()
        
        if self.backend == "sqlite":
            cursor = self.conn.cursor()
            cursor.execute(query)
            columns = [description[0] for description in cursor.description]
            data = cursor.fetchall()
            df = pd.DataFrame(data, columns=columns)
            cursor.close()
            return df
    
    
    def Identify_Twins(self, df_ID_PATIENT=None):
        '''
        Method for identifying twins.
        
        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
        
        Returns
        -------
        df_jum : DataFrame
            DataFrame containing the unique identifiers of twins population.
        '''
        
        if set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}:
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")
        
        liste_bennirpsa_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_NIR_PSA))
        
        query_JUM = f"""
           WITH MultiGem AS (
            SELECT BEN_NIR_PSA
            FROM IR_BEN_R
            WHERE BEN_NIR_PSA IN ({liste_bennirpsa_str})
            GROUP BY BEN_NIR_PSA
            HAVING COUNT(DISTINCT BEN_RNG_GEM) > 1
            )
            SELECT BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM
            FROM IR_BEN_R
            WHERE BEN_NIR_PSA IN ({liste_bennirpsa_str})
              AND BEN_NIR_PSA IN (SELECT BEN_NIR_PSA FROM MultiGem)
            ORDER BY BEN_NIR_PSA
            """
        
        df_jum = self.GetQuery(query_JUM)
       
    
        print('We have ' + str(len(np.unique(df_jum.BEN_IDT_ANO))) + ' twins in the targeted population.')

        return df_jum
    

    def loc_ccam_dcir(self, df_ID_PATIENT=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], list_CCAM=None, print_option=True):
        '''
        Method for gathering specific CCAM data from the targeted population in the DCIR.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the columns "BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM", which holds the unique identifiers of the targeted population. If None, no filter is applied on patients identifiers.
        years : list
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        list_CCAM : list
            List of CCAM codes to be retrieved. If None all CCAM codes will be returned for the targeted population. Default is None.
        print_option : bool, optional
            If True, prints the number of unique patients identified with at least one of the specified CCAM codes. Default is True.

        Returns
        -------
        df_ccam_dcir : DataFrame
            DataFrame containing CCAM codes ('CAM_PRS_IDE') from the DCIR, 
            along with their corresponding execution dates ('EXE_SOI_DTD', 'EXE_SOI_DTF'),
            for each patient ('BEN_IDT_ANO').
        '''

        if (df_ID_PATIENT is not None) and (set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}):
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        if list_CCAM is not None :
            if type(list_CCAM) != list :
                raise ValueError("list_CCAM must be a list of CCAM medical codes.")

        if (years is None) or (not isinstance(years, list)):
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")


        if isinstance(years[0], datetime):
            flxmin_year = years[0].year
            deb = years[0]
        else:
            flxmin_year = years[0]
            deb = datetime(years[0], 1, 1)

        if isinstance(years[1], datetime):
            flxmax_year = years[1].year
            end = years[1]
        else:
            flxmax_year = years[1]
            end = datetime(years[1], 12, 31)

        flxmin = datetime(flxmin_year, 1, 1)
        flxmax_date = (end + relativedelta(months=6)).replace(day=1)
        vecflx = []
        current_date = flxmin
        while current_date <= flxmax_date:
            vecflx.append(current_date.strftime("%Y-%m-%d"))
            current_date += relativedelta(months=1)

        df_ccam_dcir = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM', 'CAM_PRS_IDE', 'EXE_SOI_DTD', 'EXE_SOI_DTF'])

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


        if top_ER_PRS_F==True :

            for flux in vecflx:

                conditions = []
                if list_CCAM:
                    liste_ccam_str = ', '.join(f"'{valeur}'" for valeur in list_CCAM)
                    conditions.append(f"A.CAM_PRS_IDE IN ({liste_ccam_str})")

                if df_ID_PATIENT is not None and not df_ID_PATIENT.empty:
                    liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
                    conditions.append(f"C.BEN_IDT_ANO IN ({liste_benidtano_str})")
                    liste_bennirpsa_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_NIR_PSA))
                    conditions.append(f"C.BEN_NIR_PSA IN ({liste_bennirpsa_str})")

                conditions.append(f"B.FLX_DIS_DTD = '{flux}'")
                conditions.append(f"B.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")

                where_condition = "WHERE " + " AND ".join(conditions)

                query_CCAM_DCIR = f"""

                    SELECT  
                    C.BEN_IDT_ANO,      
                    C.BEN_NIR_PSA,
                    C.BEN_RNG_GEM,  
                    A.CAM_PRS_IDE,
                    B.EXE_SOI_DTD,
                    B.EXE_SOI_DTF
                    FROM ER_CAM_F A

                    INNER JOIN ER_PRS_F B
                    ON A.DCT_ORD_NUM = B.DCT_ORD_NUM
                    AND A.FLX_DIS_DTD = B.FLX_DIS_DTD
                    AND A.FLX_EMT_NUM = B.FLX_EMT_NUM
                    AND A.FLX_EMT_ORD = B.FLX_EMT_ORD
                    AND A.FLX_EMT_TYP = B.FLX_EMT_TYP
                    AND A.FLX_TRT_DTD = B.FLX_TRT_DTD
                    AND A.ORG_CLE_NUM = B.ORG_CLE_NUM
                    AND A.PRS_ORD_NUM = B.PRS_ORD_NUM
                    AND A.REM_TYP_AFF = B.REM_TYP_AFF

                    INNER JOIN IR_BEN_R C
                    ON B.BEN_NIR_PSA = C.BEN_NIR_PSA AND B.BEN_RNG_GEM = C.BEN_RNG_GEM

                    {where_condition}
                """
                df_flux = self.GetQuery(query_CCAM_DCIR)
                df_ccam_dcir = pd.concat([df_ccam_dcir, df_flux], ignore_index=True)
                df_ccam_dcir.drop_duplicates(inplace=True)
                df_ccam_dcir.reset_index(drop=True, inplace=True)

        else :

            conditions = []
            if list_CCAM:
                liste_ccam_str = ', '.join(f"'{valeur}'" for valeur in list_CCAM)
                conditions.append(f"A.CAM_PRS_IDE IN ({liste_ccam_str})")

            if df_ID_PATIENT is not None and not df_ID_PATIENT.empty:
                liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
                conditions.append(f"C.BEN_IDT_ANO IN ({liste_benidtano_str})")
                liste_bennirpsa_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_NIR_PSA))
                conditions.append(f"C.BEN_NIR_PSA IN ({liste_bennirpsa_str})")
                    
            conditions.append(f"B.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")

            where_condition = "WHERE " + " AND ".join(conditions)

            yy = list(range(int(flxmin_year), (int(flxmax_year)+1)))

            for year in yy :

                query_CCAM_DCIR = f"""

                    SELECT  
                    C.BEN_IDT_ANO,      
                    C.BEN_NIR_PSA,
                    C.BEN_RNG_GEM,  
                    A.CAM_PRS_IDE,
                    B.EXE_SOI_DTD,
                    B.EXE_SOI_DTF
                    FROM ER_CAM_F_{year} A

                    INNER JOIN ER_PRS_F_{year} B
                    ON A.DCT_ORD_NUM = B.DCT_ORD_NUM
                    AND A.FLX_DIS_DTD = B.FLX_DIS_DTD
                    AND A.FLX_EMT_NUM = B.FLX_EMT_NUM
                    AND A.FLX_EMT_ORD = B.FLX_EMT_ORD
                    AND A.FLX_EMT_TYP = B.FLX_EMT_TYP
                    AND A.FLX_TRT_DTD = B.FLX_TRT_DTD
                    AND A.ORG_CLE_NUM = B.ORG_CLE_NUM
                    AND A.PRS_ORD_NUM = B.PRS_ORD_NUM
                    AND A.REM_TYP_AFF = B.REM_TYP_AFF

                    INNER JOIN IR_BEN_R C
                    ON B.BEN_NIR_PSA = C.BEN_NIR_PSA AND B.BEN_RNG_GEM = C.BEN_RNG_GEM

                    {where_condition}
                """
                df_flux = self.GetQuery(query_CCAM_DCIR)
                df_ccam_dcir = pd.concat([df_ccam_dcir, df_flux], ignore_index=True)
                df_ccam_dcir.drop_duplicates(inplace=True)
                df_ccam_dcir.reset_index(drop=True, inplace=True)

        if print_option==True :
            print(str(len(np.unique(df_ccam_dcir['BEN_IDT_ANO']))) + ' patients identified using CCAM code in the DCIR.')

        return df_ccam_dcir


    def loc_ccam_pmsi(self, df_ID_PATIENT=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], list_CCAM=None, print_option=True, dev=False):
        '''
        Method for gathering specific CCAM data from the targeted population in the PMSI.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the columns "BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM", which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
        years : list
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        list_CCAM : list
            List of CCAM codes to be retrieved. If None all CCAM codes will be returned for the targeted population. Default is None.
        print_option : bool, optional
            If True, prints the number of unique patients identified with at least one of the specified CCAM codes. Default is True.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.


        Returns
        -------
        df_ccam_pmsi : DataFrame
            DataFrame containing CCAM codes ('CDC_ACT') from the PMSI, 
            along with their execution dates ('EXE_SOI_DTD', 'EXE_SOI_DTF', and 'ENT_DAT_DEL'),
            for each patient ('BEN_IDT_ANO') and specific hospital stays ('ETA_NUM', 'RSA_NUM').
        '''

        if (df_ID_PATIENT is not None) and (set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}):
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        if (years is None) or (not isinstance(years, list)): 
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")
                                
        if (list_CCAM is not None) and (type(list_CCAM) != list) :
                raise ValueError("list_CCAM must be a list of CCAM medical codes.")

        if isinstance(years[0], datetime):
            year_deb = int(str(years[0].year)[-2:])
            deb = years[0]
        else:
            year_deb = str(years[0])[-2:]
            deb = datetime(years[0], 1, 1)

        if isinstance(years[1], datetime):
            year_end = int(str(years[1].year)[-2:])
            end = years[1]
        else:
            year_end = str(years[1])[-2:]
            end = datetime(years[1], 12, 31)

        conditions = []
        if list_CCAM:
            liste_ccam_str = ', '.join(f"'{val}'" for val in list_CCAM)
            conditions.append(f"B.CDC_ACT IN ({liste_ccam_str})")
        if df_ID_PATIENT is not None and not df_ID_PATIENT.empty:
            liste_benidtano_str = ', '.join(f"'{val}'" for val in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
            conditions.append(f"C.BEN_IDT_ANO IN ({liste_benidtano_str})")
            liste_bennirpsa_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_NIR_PSA))
            conditions.append(f"C.BEN_NIR_PSA IN ({liste_bennirpsa_str})")
                
        conditions.append(f"A.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")
        where_condition = "WHERE " + " AND ".join(conditions)
    
        if dev == True :
            query_CCAM_PMSI = f"""

                SELECT DISTINCT
                    C.BEN_IDT_ANO, 
                    C.BEN_RNG_GEM, 
                    C.BEN_NIR_PSA, 
                    B.CDC_ACT, 
                    B.ENT_DAT_DEL, 
                    A.EXE_SOI_DTD, 
                    A.EXE_SOI_DTF
                FROM T_MCOaaA B

                INNER JOIN T_MCOaaC A 
                    ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                INNER JOIN IR_BEN_R C
                    ON A.NIR_ANO_17 = C.BEN_NIR_PSA
 
                {where_condition} 
            """
            df_ccam_pmsi = self.GetQuery(query_CCAM_PMSI)
            df_ccam_pmsi.drop_duplicates(inplace=True)
            df_ccam_pmsi.reset_index(drop=True, inplace=True)

        else :
            df_ccam_pmsi = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'CDC_ACT', 'ENT_DAT_DEL', 'EXE_SOI_DTD', 'EXE_SOI_DTF'])
            yy = list(range(int(year_deb), int(year_end)+1))
            for year in yy :
                table_nameA = f"T_MCO{year}A"
                table_nameC = f"T_MCO{year}C"
                query_CCAM_PMSI = f"""
                    
                    SELECT DISTINCT
                        C.BEN_IDT_ANO, 
                        C.BEN_RNG_GEM, 
                        C.BEN_NIR_PSA, 
                        B.ETA_NUM, 
                        B.RSA_NUM, 
                        B.CDC_ACT, 
                        B.ENT_DAT_DEL, 
                        A.EXE_SOI_DTD, 
                        A.EXE_SOI_DTF
                    FROM {table_nameA} B

                    INNER JOIN {table_nameC} A 
                        ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                    INNER JOIN IR_BEN_R C
                        ON A.NIR_ANO_17 = C.BEN_NIR_PSA
 
                    {where_condition} 
                """
                df_flux = self.GetQuery(query_CCAM_PMSI)
                df_ccam_pmsi = pd.concat([df_ccam_pmsi, df_flux], ignore_index=True)
                df_ccam_pmsi.drop_duplicates(inplace=True)
                df_ccam_pmsi.reset_index(drop=True, inplace=True)

        if print_option==True:
            print(str(len(np.unique(df_ccam_pmsi['BEN_IDT_ANO']))) + ' patient identified using CCAM code in the PMSI.')
        
        return df_ccam_pmsi
    
    
    
    def loc_icd10_pmsi(self, df_ID_PATIENT=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], list_ICD10=None, print_option=True, dev=False):
        '''
        Method for gathering specific ICD-10 data from the targeted population in the PMSI.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the columns "BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM", which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
        years : list
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for ICD-10 codes. By default 1st of January 2020 and 31 of December 2020.
        list_ICD10 : list
            List of ICD-10 codes to be retrieved. If None all ICD-10 codes will be returned for the targeted population. Default is None.
        print_option : bool, optional
            If True, prints the number of unique patients identified with at least one of the specified ICD-10 codes. Default is True.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.

        Returns
        -------
        df_icd10_pmsi : DataFrame
            DataFrame containing ICD-10 codes ('DGN_PAL', 'DGN_REL', 'ASS_DGN') from the PMSI, 
            along with their execution dates ('EXE_SOI_DTD', 'EXE_SOI_DTF'), 
            for each patient ('BEN_IDT_ANO') and specific hospital stays ('ETA_NUM', 'RSA_NUM').
        '''

        if (df_ID_PATIENT is not None) and (set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}):
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        if (years is None) or (not isinstance(years, list)): 
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")    

        if (list_ICD10 is not None) and (type(list_ICD10) != list) :
                raise ValueError("list_ICD10 must be a list of ICD-10 medical codes.")

    
        if isinstance(years[0], datetime):
            year_deb = int(str(years[0].year)[-2:])
            deb = years[0]
        else:
            year_deb = str(years[0])[-2:]
            deb = datetime(years[0], 1, 1)

        if isinstance(years[1], datetime):
            year_end = int(str(years[1].year)[-2:])
            end = years[1]
        else:
            year_end = str(years[1])[-2:]
            end = datetime(years[1], 12, 31)


        condition_icd10_1, condition_icd10_2 = [], []

        if df_ID_PATIENT is not None and not df_ID_PATIENT.empty:
            liste_benidtano_str = ', '.join(f"'{val}'" for val in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
            condition_icd10_1.append(f"C.BEN_IDT_ANO IN ({liste_benidtano_str})")
            condition_icd10_2.append(f"C.BEN_IDT_ANO IN ({liste_benidtano_str})")
            
            liste_bennirpsa_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_NIR_PSA))
            condition_icd10_1.append(f"C.BEN_NIR_PSA IN ({liste_bennirpsa_str})")
            condition_icd10_2.append(f"C.BEN_NIR_PSA IN ({liste_bennirpsa_str})")

        if list_ICD10:
            liste_icd10_str = ', '.join(f"'{valeur}'" for valeur in list_ICD10)
            condition_icd10_1.append(f"(B.DGN_PAL IN ({liste_icd10_str}) OR B.DGN_REL IN ({liste_icd10_str}))")
            condition_icd10_2.append(f"E.ASS_DGN IN ({liste_icd10_str})")
       
        condition_icd10_1.append(f"A.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")
        condition_icd10_2.append(f"A.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")

        where_condition_1 = "WHERE " + " AND ".join(condition_icd10_1)
        where_condition_2 = "WHERE " + " AND ".join(condition_icd10_2)
    
        if dev == True :
            query_ICD10_PMSI = f"""
    
                SELECT DISTINCT
                C.BEN_IDT_ANO, 
                C.BEN_RNG_GEM, 
                C.BEN_NIR_PSA, 
                B.DGN_PAL, 
                B.DGN_REL, 
                NULL AS ASS_DGN,
                A.EXE_SOI_DTD, 
                A.EXE_SOI_DTF
                FROM T_MCOaaB B

                INNER JOIN T_MCOaaC A 
                ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                INNER JOIN IR_BEN_R C
                ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                {where_condition_1} 
                
                UNION ALL
                
                SELECT 
                C.BEN_IDT_ANO,
                C.BEN_RNG_GEM, 
                C.BEN_NIR_PSA, 
                NULL AS DGN_PAL,
                NULL AS DGN_REL,
                E.ASS_DGN, 
                A.EXE_SOI_DTD, 
                A.EXE_SOI_DTF
                FROM T_MCOaaD E

                INNER JOIN T_MCOaaC A 
                ON E.ETA_NUM = A.ETA_NUM AND E.RSA_NUM = A.RSA_NUM

                INNER JOIN IR_BEN_R C
                ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                {where_condition_2}
            """
            df_icd10_pmsi = self.GetQuery(query_ICD10_PMSI) 
            df_icd10_pmsi.drop_duplicates(inplace=True)
            df_icd10_pmsi.reset_index(drop=True, inplace=True)


        else :
            df_icd10_pmsi = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'DGN_PAL', 'DGN_REL', 'ASS_DGN', 'EXE_SOI_DTD', 'EXE_SOI_DTF'])
            yy = list(range(int(year_deb), int(year_end)+1))
            for year in yy :
                table_nameB = f"T_MCO{year}B"
                table_nameC = f"T_MCO{year}C"
                table_nameD = f"T_MCO{year}D"
                query_ICD10_PMSI = f"""
                    
                    SELECT DISTINCT
                    C.BEN_IDT_ANO, 
                    C.BEN_RNG_GEM, 
                    C.BEN_NIR_PSA, 
                    B.DGN_PAL, 
                    B.DGN_REL, 
                    NULL AS ASS_DGN,
                    A.EXE_SOI_DTD, 
                    A.EXE_SOI_DTF
                    FROM {table_nameB} B

                    INNER JOIN {table_nameC} A 
                    ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                    INNER JOIN IR_BEN_R C
                    ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                    {where_condition_1}
                    
                    UNION ALL
                    
                    SELECT 
                    C.BEN_IDT_ANO,
                    C.BEN_RNG_GEM, 
                    C.BEN_NIR_PSA, 
                    NULL AS DGN_PAL,
                    NULL AS DGN_REL,
                    E.ASS_DGN, 
                    A.EXE_SOI_DTD, 
                    A.EXE_SOI_DTF
                    FROM {table_nameD} E

                    INNER JOIN {table_nameC} A 
                    ON E.ETA_NUM = A.ETA_NUM AND E.RSA_NUM = A.RSA_NUM

                    INNER JOIN IR_BEN_R C
                    ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                    {where_condition_2}
                """
                df_flux = self.GetQuery(query_ICD10_PMSI) 
                df_icd10_pmsi = pd.concat([df_icd10_pmsi, df_flux], ignore_index=True)
                df_icd10_pmsi.drop_duplicates(inplace=True)
                df_icd10_pmsi.reset_index(drop=True, inplace=True)

        if print_option==True :
            print(str(len(np.unique(df_icd10_pmsi['BEN_IDT_ANO']))) + ' patients identified using ICD10 code in the PMSI.')
        
        return df_icd10_pmsi
    

    def loc_ucd_pmsi(self, df_ID_PATIENT=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], list_UCD=None, print_option=True, dev=False):
        '''
        Method for gathering specific UCD data from the targeted population in the PMSI.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
        DataFrame containing the columns "BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM", which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
        years : list
        List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for UCD codes. By default 1st of January 2020 and 31 of December 2020.
        list_UCD : list
        List of UCD codes to be retrieved. If None all UCD codes will be returned for the targeted population. Default is None.
        print_option : bool, optional
        If True, prints the number of unique patients identified with at least one of the specified UCD codes. Default is True.
        dev : bool, optional
        If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.

        Returns
        -------
        df_ucd_pmsi : DataFrame
        DataFrame containing UCD codes ('UCD_UCD_COD', 'UCD_COD') from the PMSI, 
        along with their execution dates ('EXE_SOI_DTD', 'EXE_SOI_DTF' and 'DELAI'), 
        for each patient ('BEN_IDT_ANO') and specific hospital stays ('ETA_NUM', 'RSA_NUM').
        '''

        if (df_ID_PATIENT is not None) and (set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}):
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        if list_UCD is not None :
            if type(list_UCD) != list :
                raise ValueError("list_UCD must be a list of UCD medical codes.")

        if (years is None) or (not isinstance(years, list)): 
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")

        if isinstance(years[0], datetime):
            year_deb = int(str(years[0].year)[-2:])
            deb = years[0]
        else:
            year_deb = str(years[0])[-2:]
            deb = datetime(years[0], 1, 1)

        if isinstance(years[1], datetime):
            year_end = int(str(years[1].year)[-2:])
            end = years[1]
        else:
            year_end = str(years[1])[-2:]
            end = datetime(years[1], 12, 31)

        conditions = []
        if list_UCD:
            liste_ucd_str = ', '.join(f"'{val}'" for val in list_UCD)
            conditions.append(f"RIGHT(B.UCD_UCD_COD, 7) IN ({liste_ucd_str})")
        if df_ID_PATIENT is not None and not df_ID_PATIENT.empty:
            liste_benidtano_str = ', '.join(f"'{val}'" for val in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
            conditions.append(f"C.BEN_IDT_ANO IN ({liste_benidtano_str})")
            liste_bennirpsa_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_NIR_PSA))
            conditions.append(f"C.BEN_NIR_PSA IN ({liste_bennirpsa_str})")
        conditions.append(f"A.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")
        where_condition = "WHERE " + " AND ".join(conditions)        

        if dev == True :
            query_UCD_PMSIMED = f"""

                SELECT DISTINCT
                C.BEN_IDT_ANO, 
                C.BEN_RNG_GEM, 
                C.BEN_NIR_PSA, 
                B.UCD_UCD_COD,
                RIGHT(B.UCD_UCD_COD, 7) AS COD_UCD,
                D.PHA_ATC_CLA,
                D.PHA_ATC_LIB,
                D.PHA_ATC_C03,
                D.PHA_ATC_L03,
                A.EXE_SOI_DTD, 
                A.EXE_SOI_DTF
                FROM T_MCOaaMED B

                INNER JOIN T_MCOaaC A 
                ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                INNER JOIN IR_BEN_R C
                ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                LEFT JOIN IR_PHA_R D
                ON RIGHT(B.UCD_UCD_COD, 7) = D.PHA_CIP_UCD

                {where_condition}
            """
            df_ucd_pmsi_MED = self.GetQuery(query_UCD_PMSIMED)
            df_ucd_pmsi_MED.drop_duplicates(inplace=True)
            df_ucd_pmsi_MED.reset_index(drop=True, inplace=True)
                
            query_UCD_PMSIFH = f"""

                SELECT DISTINCT
                C.BEN_IDT_ANO, 
                C.BEN_RNG_GEM, 
                C.BEN_NIR_PSA, 
                B.UCD_UCD_COD,
                RIGHT(B.UCD_UCD_COD, 7) AS COD_UCD,
                D.PHA_ATC_CLA,
                D.PHA_ATC_LIB,
                D.PHA_ATC_C03,
                D.PHA_ATC_L03,
                A.EXE_SOI_DTD, 
                A.EXE_SOI_DTF
                FROM T_MCOaaFH B

                INNER JOIN T_MCOaaC A 
                ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                INNER JOIN IR_BEN_R C
                ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                LEFT JOIN IR_PHA_R D
                ON RIGHT(B.UCD_UCD_COD, 7) = D.PHA_CIP_UCD

                {where_condition}
            """
            df_ucd_pmsi_FH = self.GetQuery(query_UCD_PMSIFH)
            df_ucd_pmsi = pd.concat([df_ucd_pmsi_MED, df_ucd_pmsi_FH], ignore_index=True)
            df_ucd_pmsi.drop_duplicates(inplace=True)
            df_ucd_pmsi.reset_index(drop=True, inplace=True)


        else :

            df_ucd_pmsi = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'UCD_UCD_COD', 'COD_UCD', 'PHA_ATC_CLA', 'PHA_ATC_LIB', 'PHA_ATC_C03', 'PHA_ATC_L03', 'EXE_SOI_DTD', 'EXE_SOI_DTF'])
            yy = list(range(int(year_deb), int(year_end)+1))

            for year in yy :

                table_nameMED = f"T_MCO{year}MED"
                table_nameFH = f"T_MCO{year}FH"
                table_nameC = f"T_MCO{year}C"

                query_UCD_PMSIMED = f"""

                    SELECT DISTINCT
                    C.BEN_IDT_ANO, 
                    C.BEN_RNG_GEM, 
                    C.BEN_NIR_PSA, 
                    B.UCD_UCD_COD,
                    RIGHT(B.UCD_UCD_COD, 7) AS COD_UCD,
                    D.PHA_ATC_CLA,
                    D.PHA_ATC_LIB,
                    D.PHA_ATC_C03,
                    D.PHA_ATC_L03,
                    A.EXE_SOI_DTD, 
                    A.EXE_SOI_DTF
                    FROM {table_nameMED} B

                    INNER JOIN {table_nameC} A 
                    ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                    INNER JOIN IR_BEN_R C
                    ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                    LEFT JOIN IR_PHA_R D
                    ON RIGHT(B.UCD_UCD_COD, 7) = D.PHA_CIP_UCD

                    {where_condition}
                """
                df_flux_MED = self.GetQuery(query_UCD_PMSIMED)
                df_ucd_pmsi = pd.concat([df_ucd_pmsi, df_flux_MED], ignore_index=True)
                df_ucd_pmsi.drop_duplicates(inplace=True)
                df_ucd_pmsi.reset_index(drop=True, inplace=True)

                query_UCD_PMSIFH = f"""

                    SELECT DISTINCT
                    C.BEN_IDT_ANO, 
                    C.BEN_RNG_GEM, 
                    C.BEN_NIR_PSA, 
                    B.UCD_UCD_COD,
                    RIGHT(B.UCD_UCD_COD, 7) AS COD_UCD,
                    D.PHA_ATC_CLA,
                    D.PHA_ATC_LIB,
                    D.PHA_ATC_C03,
                    D.PHA_ATC_L03,
                    A.EXE_SOI_DTD, 
                    A.EXE_SOI_DTF
                    FROM {table_nameFH} B

                    INNER JOIN {table_nameC} A 
                    ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                    INNER JOIN IR_BEN_R C
                    ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                    LEFT JOIN IR_PHA_R D
                    ON RIGHT(B.UCD_UCD_COD, 7) = D.PHA_CIP_UCD

                    {where_condition}
                """
                df_flux_FH = self.GetQuery(query_UCD_PMSIFH)
                df_ucd_pmsi = pd.concat([df_ucd_pmsi, df_flux_FH], ignore_index=True)
                df_ucd_pmsi.drop_duplicates(inplace=True)
                df_ucd_pmsi.reset_index(drop=True, inplace=True)

        if print_option==True :
            print(str(len(np.unique(df_ucd_pmsi['BEN_IDT_ANO']))) + ' patients identified using UCD code in the PMSI.')

        return df_ucd_pmsi


    def loc_cip_dcir(self, df_ID_PATIENT=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], list_CIP13=None, print_option=True):
        '''
        Method for gathering specific UCD data from the targeted population in the DCIR.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the columns "BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM", which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
        years : list
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CIP codes. By default 1st of January 2020 and 31 of December 2020.
        list_CIP13 : list
            List of CIP-13 codes to be retrieved. If None all CIP-13 codes will be returned for the targeted population. Default is None.
        print_option : bool, optional
            If True, prints the number of unique patients identified with at least one of the specified CIP-13 codes. Default is True.
        
        Returns
        -------
        df_cip_dcir : DataFrame
            DataFrame containing CIP-13 codes ('PHA_CIP_C13') from the DCIR, 
            with their equivalence in ATC coding ('PHA_ATC_CLA', 'PHA_ATC_LIB'), 
            along with their corresponding execution dates ('EXE_SOI_DTD', 'EXE_SOI_DTF'), 
            for each patient ('BEN_IDT_ANO').
        '''
        
        if (df_ID_PATIENT is not None) and (set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}):
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        if list_CIP13 is not None :
            if type(list_CIP13) != list :
                raise ValueError("list_CIP13 must be a list of CIP-13 medical codes.")

        if (years is None) or (not isinstance(years, list)): 
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")

        if isinstance(years[0], datetime):
            flxmin_year = years[0].year
            deb = years[0]
        else:
            flxmin_year = years[0]
            deb = datetime(years[0], 1, 1)

        if isinstance(years[1], datetime):
            flxmax_year = years[1].year
            end = years[1]
        else:
            flxmax_year = years[1]
            end = datetime(years[1], 12, 31)

        flxmin = datetime(flxmin_year, 1, 1)
        flxmax_date = (end + relativedelta(months=6)).replace(day=1)
        vecflx = []
        current_date = flxmin
        while current_date <= flxmax_date:
            vecflx.append(current_date.strftime("%Y-%m-%d"))
            current_date += relativedelta(months=1)

        df_cip_dcir = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM', 'EXE_SOI_DTD', 'EXE_SOI_DTF', 'PHA_CIP_C13', 'PHA_ATC_CLA', 'PHA_ATC_LIB', 'PHA_ATC_C03', 'PHA_ATC_L03'])


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

            
        if top_ER_PRS_F==True :

            for flux in vecflx:

                conditions = []

                if list_CIP13:
                    liste_cip13_str = ', '.join(f"'{valeur}'" for valeur in list_CIP13)
                    conditions.append(f"D.PHA_PRS_C13 IN ({liste_cip13_str})")

                if df_ID_PATIENT is not None and not df_ID_PATIENT.empty:
                    liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
                    conditions.append(f"B.BEN_IDT_ANO IN ({liste_benidtano_str})")
                    liste_bennirpsa_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_NIR_PSA))
                    conditions.append(f"B.BEN_NIR_PSA IN ({liste_bennirpsa_str})")

                conditions.append(f"A.FLX_DIS_DTD = '{flux}'")
                conditions.append(f"A.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")

                where_condition = "WHERE " + " AND ".join(conditions)

                query_CIP13_DCIR = f"""

                    SELECT DISTINCT
                    B.BEN_IDT_ANO,
                    A.BEN_NIR_PSA, 
                    A.BEN_RNG_GEM, 
                    A.EXE_SOI_DTD,
                    A.EXE_SOI_DTF,
                    D.PHA_PRS_C13 AS PHA_CIP_C13,
                    E.PHA_ATC_CLA,
                    E.PHA_ATC_LIB,
                    E.PHA_ATC_C03,
                    E.PHA_ATC_L03
                    FROM ER_PRS_F A

                    INNER JOIN IR_BEN_R B
                    ON (A.BEN_NIR_PSA = B.BEN_NIR_PSA AND A.BEN_RNG_GEM = B.BEN_RNG_GEM)

                    INNER JOIN ER_PHA_F D ON (
                    A.FLX_DIS_DTD = D.FLX_DIS_DTD 
                    AND A.FLX_TRT_DTD = D.FLX_TRT_DTD
                    AND A.FLX_EMT_TYP = D.FLX_EMT_TYP
                    AND A.FLX_EMT_NUM = D.FLX_EMT_NUM
                    AND A.FLX_EMT_ORD = D.FLX_EMT_ORD
                    AND A.ORG_CLE_NUM = D.ORG_CLE_NUM
                    AND A.DCT_ORD_NUM = D.DCT_ORD_NUM
                    AND A.PRS_ORD_NUM = D.PRS_ORD_NUM
                    AND A.REM_TYP_AFF = D.REM_TYP_AFF
                    )

                    LEFT JOIN IR_PHA_R E ON 
                    D.PHA_PRS_C13 = E.PHA_CIP_C13

                    {where_condition}
                """
                df_flux = self.GetQuery(query_CIP13_DCIR)
                df_cip_dcir = pd.concat([df_cip_dcir, df_flux], ignore_index=True)
                df_cip_dcir.drop_duplicates(inplace=True)
                df_cip_dcir.reset_index(drop=True, inplace=True)


        else :

            conditions = []

            if list_CIP13:
                liste_cip13_str = ', '.join(f"'{valeur}'" for valeur in list_CIP13)
                conditions.append(f"D.PHA_PRS_C13 IN ({liste_cip13_str})")

            if df_ID_PATIENT is not None and not df_ID_PATIENT.empty:
                liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
                conditions.append(f"B.BEN_IDT_ANO IN ({liste_benidtano_str})")
                liste_bennirpsa_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_NIR_PSA))
                conditions.append(f"B.BEN_NIR_PSA IN ({liste_bennirpsa_str})")

            conditions.append(f"A.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")

            where_condition = "WHERE " + " AND ".join(conditions)

            yy = list(range(int(flxmin_year), (int(flxmax_year)+1)))

            for year in yy :

                query_CIP13_DCIR = f"""

                    SELECT DISTINCT
                    B.BEN_IDT_ANO,
                    A.BEN_NIR_PSA, 
                    A.BEN_RNG_GEM, 
                    A.EXE_SOI_DTD,
                    A.EXE_SOI_DTF,
                    D.PHA_PRS_C13 AS PHA_CIP_C13,
                    E.PHA_ATC_CLA,
                    E.PHA_ATC_LIB,
                    E.PHA_ATC_C03,
                    E.PHA_ATC_L03
                    FROM ER_PRS_F_{year} A

                    INNER JOIN IR_BEN_R B
                    ON (A.BEN_NIR_PSA = B.BEN_NIR_PSA AND A.BEN_RNG_GEM = B.BEN_RNG_GEM)

                    INNER JOIN ER_PHA_F_{year} D ON (
                    A.FLX_DIS_DTD = D.FLX_DIS_DTD 
                    AND A.FLX_TRT_DTD = D.FLX_TRT_DTD
                    AND A.FLX_EMT_TYP = D.FLX_EMT_TYP
                    AND A.FLX_EMT_NUM = D.FLX_EMT_NUM
                    AND A.FLX_EMT_ORD = D.FLX_EMT_ORD
                    AND A.ORG_CLE_NUM = D.ORG_CLE_NUM
                    AND A.DCT_ORD_NUM = D.DCT_ORD_NUM
                    AND A.PRS_ORD_NUM = D.PRS_ORD_NUM
                    AND A.REM_TYP_AFF = D.REM_TYP_AFF
                    )

                    LEFT JOIN IR_PHA_R E ON 
                    D.PHA_PRS_C13 = E.PHA_CIP_C13

                    {where_condition}
                """
                df_flux = self.GetQuery(query_CIP13_DCIR)
                df_cip_dcir = pd.concat([df_cip_dcir, df_flux], ignore_index=True)
                df_cip_dcir.drop_duplicates(inplace=True)
                df_cip_dcir.reset_index(drop=True, inplace=True)


        if print_option==True :
            print(str(len(np.unique(df_cip_dcir['BEN_IDT_ANO']))) + ' patients identified using CIP code in the DCIR.')
    
        return df_cip_dcir
    
    
    def loc_ucd_dcir(self, df_ID_PATIENT=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], list_UCD=None, print_option=True):
        '''
        Method for gathering specific UCD data from the targeted population in the DCIR.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
        DataFrame containing the columns "BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM", which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
        years : list
        List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CIP codes. By default 1st of January 2020 and 31 of December 2020.
        list_UCD : list
        List of UCD codes to be retrieved. If None all UCD codes will be returned for the targeted population. Default is None.
        print_option : bool, optional
        If True, prints the number of unique patients identified with at least one of the specified CIP-13 codes. Default is True.

        Returns
        -------
        df_ucd_dcir : DataFrame
        DataFrame containing UCD codes from the DCIR, 
        with their equivalence in ATC coding ('PHA_ATC_CLA', 'PHA_ATC_LIB'), 
        along with their corresponding execution dates ('EXE_SOI_DTD', 'EXE_SOI_DTF'), 
        for each patient ('BEN_IDT_ANO').
        '''
        
        if (df_ID_PATIENT is not None) and (set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}):
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        if list_UCD is not None :
            if type(list_UCD) != list :
                raise ValueError("list_UCD must be a list of UCD medical codes.")

        if (years is None) or (not isinstance(years, list)): 
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")

        if isinstance(years[0], datetime):
            flxmin_year = years[0].year
            deb = years[0]
        else:
            flxmin_year = years[0]
            deb = datetime(years[0], 1, 1)

        if isinstance(years[1], datetime):
            flxmax_year = years[1].year
            end = years[1]
        else:
            flxmax_year = years[1]
            end = datetime(years[1], 12, 31)

        flxmin = datetime(flxmin_year, 1, 1)
        flxmax_date = (end + relativedelta(months=6)).replace(day=1)
        vecflx = []
        current_date = flxmin
        while current_date <= flxmax_date:
            vecflx.append(current_date.strftime("%Y-%m-%d"))
            current_date += relativedelta(months=1)

        df_ucd_dcir = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM', 'EXE_SOI_DTD', 'EXE_SOI_DTF', 'UCD_UCD_COD', 'COD_UCD', 'PHA_ATC_CLA', 'PHA_ATC_LIB','PHA_ATC_C03', 'PHA_ATC_L03'])


        # top_ER_PRS_F = True
        # if self.backend == 'sqlite':
        #     cursor = self.conn.cursor()
        #     cursor.execute(
        #     "SELECT name FROM sqlite_master WHERE type='table' AND name=?;",
        #     ('E_PRS_F',)
        #     )
        #     top_ER_PRS_F = cursor.fetchone() is not None


        if self.backend == 'spark':
            top_ER_PRS_F = self.conn.catalog.tableExists('ER_PRS_F')


        if top_ER_PRS_F==True :

            for flux in vecflx:

                conditions = []
                if list_UCD:
                    liste_ucd_str = ', '.join(f"'{val}'" for val in list_UCD)
                    conditions.append(f"RIGHT(B.UCD_UCD_COD, 7) IN ({liste_ucd_str})")
                if df_ID_PATIENT is not None and not df_ID_PATIENT.empty:
                    liste_benidtano_str = ', '.join(f"'{val}'" for val in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
                    conditions.append(f"C.BEN_IDT_ANO IN ({liste_benidtano_str})")
                    liste_bennirpsa_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_NIR_PSA))
                    conditions.append(f"C.BEN_NIR_PSA IN ({liste_bennirpsa_str})")
                    
                conditions.append(f"A.FLX_DIS_DTD = '{flux}'")
                conditions.append(f"A.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")

                where_condition = "WHERE " + " AND ".join(conditions)        

                query_UCD_ER_UCD_F = f"""

                    SELECT DISTINCT
                        C.BEN_IDT_ANO, 
                        C.BEN_RNG_GEM, 
                        C.BEN_NIR_PSA, 
                        B.UCD_UCD_COD,
                        RIGHT(B.UCD_UCD_COD, 7) AS COD_UCD,
                        D.PHA_ATC_CLA,
                        D.PHA_ATC_LIB,
                        D.PHA_ATC_C03,
                        D.PHA_ATC_L03,
                        A.EXE_SOI_DTD, 
                        A.EXE_SOI_DTF
                    FROM ER_UCD_F B

                    INNER JOIN ER_PRS_F A ON (
                            A.FLX_DIS_DTD = B.FLX_DIS_DTD 
                            AND A.FLX_TRT_DTD = B.FLX_TRT_DTD
                            AND A.FLX_EMT_TYP = B.FLX_EMT_TYP
                            AND A.FLX_EMT_NUM = B.FLX_EMT_NUM
                            AND A.FLX_EMT_ORD = B.FLX_EMT_ORD
                            AND A.ORG_CLE_NUM = B.ORG_CLE_NUM
                            AND A.DCT_ORD_NUM = B.DCT_ORD_NUM
                            AND A.PRS_ORD_NUM = B.PRS_ORD_NUM
                            AND A.REM_TYP_AFF = B.REM_TYP_AFF
                        )

                    INNER JOIN IR_BEN_R C
                        ON A.BEN_NIR_PSA = C.BEN_NIR_PSA

                    LEFT JOIN IR_PHA_R D
                        ON RIGHT(B.UCD_UCD_COD, 7) = D.PHA_CIP_UCD

                    {where_condition}

                """

                df_flux = self.GetQuery(query_UCD_ER_UCD_F)
                df_ucd_dcir = pd.concat([df_ucd_dcir, df_flux], ignore_index=True)
                df_ucd_dcir.drop_duplicates(inplace=True)
                df_ucd_dcir.reset_index(drop=True, inplace=True)


        else :

            conditions = []
            if list_UCD:
                liste_ucd_str = ', '.join(f"'{val}'" for val in list_UCD)
                conditions.append(f"RIGHT(B.UCD_UCD_COD, 7) IN ({liste_ucd_str})")
            if df_ID_PATIENT is not None and not df_ID_PATIENT.empty:
                liste_benidtano_str = ', '.join(f"'{val}'" for val in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
                conditions.append(f"C.BEN_IDT_ANO IN ({liste_benidtano_str})")
                liste_bennirpsa_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_NIR_PSA))
                conditions.append(f"C.BEN_NIR_PSA IN ({liste_bennirpsa_str})")

            conditions.append(f"A.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")

            where_condition = "WHERE " + " AND ".join(conditions)      

            yy = list(range(int(flxmin_year), (int(flxmax_year)+1)))
            
            for year in yy :

                query_UCD_ER_UCD_F = f"""

                    SELECT DISTINCT
                    C.BEN_IDT_ANO, 
                    C.BEN_RNG_GEM, 
                    C.BEN_NIR_PSA, 
                    B.UCD_UCD_COD,
                    RIGHT(B.UCD_UCD_COD, 7) AS COD_UCD,
                    D.PHA_ATC_CLA,
                    D.PHA_ATC_LIB,
                    D.PHA_ATC_C03,
                    D.PHA_ATC_L03,
                    A.EXE_SOI_DTD, 
                    A.EXE_SOI_DTF
                    FROM ER_UCD_F_{year} B

                    INNER JOIN ER_PRS_F_{year} A ON (
                    A.FLX_DIS_DTD = B.FLX_DIS_DTD 
                    AND A.FLX_TRT_DTD = B.FLX_TRT_DTD
                    AND A.FLX_EMT_TYP = B.FLX_EMT_TYP
                    AND A.FLX_EMT_NUM = B.FLX_EMT_NUM
                    AND A.FLX_EMT_ORD = B.FLX_EMT_ORD
                    AND A.ORG_CLE_NUM = B.ORG_CLE_NUM
                    AND A.DCT_ORD_NUM = B.DCT_ORD_NUM
                    AND A.PRS_ORD_NUM = B.PRS_ORD_NUM
                    AND A.REM_TYP_AFF = B.REM_TYP_AFF
                    )

                    INNER JOIN IR_BEN_R C
                    ON A.BEN_NIR_PSA = C.BEN_NIR_PSA

                    LEFT JOIN IR_PHA_R D
                    ON RIGHT(B.UCD_UCD_COD, 7) = D.PHA_CIP_UCD

                    {where_condition}

                """

                df_flux = self.GetQuery(query_UCD_ER_UCD_F)
                df_ucd_dcir = pd.concat([df_ucd_dcir, df_flux], ignore_index=True)
                df_ucd_dcir.drop_duplicates(inplace=True)
                df_ucd_dcir.reset_index(drop=True, inplace=True)

        if print_option==True :
            print(str(len(np.unique(df_ucd_dcir['BEN_IDT_ANO']))) + ' patients identified using UCD code in the DCIR.')

        return df_ucd_dcir

    


    def loc_atc_pmsi(self, list_ATC, df_ID_PATIENT=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], print_option=True, dev=False):
        '''
        Method for gathering specific ATC data from the targeted population in the PMSI.

        Parameters
        ----------
        list_ATC : list
        List of ATC codes to be retrieved. If None all CCAM codes will be returned for the targeted population.
        df_ID_PATIENT : DataFrame
        DataFrame containing the columns "BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM", which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
        years : list
        List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        print_option : bool, optional
        If True, prints the number of unique patients identified with at least one of the specified CCAM codes. Default is True.
        dev : bool, optional
        If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.


        Returns
        -------
        df_atc_pmsi : DataFrame
        DataFrame containing ATC codes and wording ('PHA_ATC_CLA', 'PHA_ATC_LIB') from the PMSI,
        with their correspondance in UCD ('PHA_CIP_UCD') 
        along with their execution dates ('EXE_SOI_DTD', 'EXE_SOI_DTF'),
        for each patient ('BEN_IDT_ANO') and specific hospital stays ('ETA_NUM', 'RSA_NUM').
        '''

        if (df_ID_PATIENT is not None) and (set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}):
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        if (list_ATC is None) or (type(list_ATC) != list):
            raise ValueError("list_ATC must be a list of ATC medical codes.")

        if (years is None) or (not isinstance(years, list)): 
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")

        if isinstance(years[0], datetime):
            year_deb = int(str(years[0].year)[-2:])
            deb = years[0]
        else:
            year_deb = str(years[0])[-2:]
            deb = datetime(years[0], 1, 1)

        if isinstance(years[1], datetime):
            year_end = int(str(years[1].year)[-2:])
            end = years[1]
        else:
            year_end = str(years[1])[-2:]
            end = datetime(years[1], 12, 31)

        liste_atc_str = ', '.join(f"'{valeur}'" for valeur in list_ATC)

        query_ATC = f"""

            SELECT DISTINCT
            PHA_ATC_CLA,
            PHA_ATC_LIB,
            PHA_CIP_C13,
            PHA_CIP_UCD
            FROM IR_PHA_R

            WHERE PHA_ATC_CLA IN ({liste_atc_str})
        """ 
        df_atc_cip_ucd = self.GetQuery(query_ATC)

        list_UCD = df_atc_cip_ucd.PHA_CIP_UCD.tolist()
        df_atc_pmsi = self.loc_ucd_pmsi(df_ID_PATIENT, years=years, list_UCD=list_UCD, print_option=print_option, dev=dev)
        df_atc_pmsi.drop_duplicates(inplace=True)
        df_atc_pmsi.reset_index(drop=True, inplace=True)

        return df_atc_pmsi


    def loc_atc_dcir(self, list_ATC, df_ID_PATIENT=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], print_option=True):
        '''
        Method for gathering specific ATC data from the targeted population in the DCIR.

        Parameters
        ----------
        list_ATC : list
        List of ATC codes to be retrieved. If None all CCAM codes will be returned for the targeted population.
        df_ID_PATIENT : DataFrame
        DataFrame containing the columns "BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM", which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
        years : list
        List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        print_option : bool, optional
        If True, prints the number of unique patients identified with at least one of the specified CCAM codes. Default is True.

        Returns
        -------
        df_atc_pmsi : DataFrame
        DataFrame containing ATC codes and wording ('PHA_ATC_CLA', 'PHA_ATC_LIB') from the DCIR,
        with their correspondance in UCD ('PHA_CIP_UCD') and CIP13 ('PHA_CIP_C13')
        along with their execution dates ('EXE_SOI_DTD', 'EXE_SOI_DTF'),
        for each patient ('BEN_IDT_ANO').
        '''

        if (df_ID_PATIENT is not None) and (set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}):
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        if (list_ATC is None) or (type(list_ATC) != list):
            raise ValueError("list_ATC must be a list of ATC medical codes.")

        if (years is None) or (not isinstance(years, list)): 
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")

        if isinstance(years[0], datetime):
            year_deb = int(str(years[0].year)[-2:])
            deb = years[0]
        else:
            year_deb = str(years[0])[-2:]
            deb = datetime(years[0], 1, 1)

        if isinstance(years[1], datetime):
            year_end = int(str(years[1].year)[-2:])
            end = years[1]
        else:
            year_end = str(years[1])[-2:]
            end = datetime(years[1], 12, 31)

        liste_atc_str = ', '.join(f"'{valeur}'" for valeur in list_ATC)

        query_ATC = f"""

            SELECT DISTINCT
            PHA_ATC_CLA,
            PHA_ATC_LIB,
            PHA_CIP_C13,
            PHA_CIP_UCD
            FROM IR_PHA_R

            WHERE PHA_ATC_CLA IN ({liste_atc_str})
        """ 
        df_atc_cip_ucd = self.GetQuery(query_ATC)

        list_UCD = df_atc_cip_ucd.PHA_CIP_UCD.tolist()
        list_CIP13 = df_atc_cip_ucd.PHA_CIP_C13.tolist()
        df_ucd_dcir = self.loc_ucd_dcir(df_ID_PATIENT, years=years, list_UCD=list_UCD, print_option=print_option)
        df_cip_dcir = self.loc_cip_dcir(df_ID_PATIENT, years=years, list_CIP13=list_CIP13, print_option=print_option)

        df_atc_dcir = pd.concat([df_ucd_dcir, df_cip_dcir]).drop_duplicates().reset_index(drop=True)
        df_atc_dcir.drop_duplicates(inplace=True)
        df_atc_dcir.reset_index(drop=True, inplace=True)
        
        return df_atc_dcir




    def Get_records(self, df_ID_PATIENT, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], list_CCAM=None, list_ICD10=None, list_UCD=None, list_CIP13=None, list_ATC=None, export=False, path='', dev=False):
        '''
        Method for getting patients records in a DataFrame and optionaly exporting it in pickle format.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the columns "BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM", which holds the unique identifiers of the targeted population.
        years : list
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for UCD codes. By default 1st of January 2020 and 31 of December 2020.
        list_CCAM : list
            List of CCAM codes to be retrieved. If None all CCAM codes will be returned for the targeted population. Default is None.
        list_ICD10 : list
            List of ICD-10 codes to be retrieved. If None all ICD-10 codes will be returned for the targeted population. Default is None.
        list_UCD : list
            List of UCD codes to be retrieved. If None all UCD codes will be returned for the targeted population. Default is None.
        list_CIP13 : list
            List of CIP-13 codes to be retrieved. If None all CIP-13 codes will be returned for the targeted population. Default is None.
        list_ATC : list
            List of ATC codes to be retrieved. If None ATC  codes will NOT be returned for the targeted population. Default is None.
        export : bool, optional
            If True, the resulting dataframe of records will be exported in pickle format. Default is False.
        path : str, optional
            The path where to export the dataframe.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.

        Returns
        -------
        df_records : DataFrame
            Medical records of patients.
        '''

        if (df_ID_PATIENT is not None) and (set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}):
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        if (list_CCAM is not None) and (type(list_CCAM) != list) :
            raise ValueError("list_CCAM must be a list of CCAM medical codes.")

        if (list_ICD10 is not None) and (type(list_ICD10) != list) :
            raise ValueError("list_ICD10 must be a list of ICD-10 medical codes.")

        if (list_UCD is not None) and (type(list_UCD) != list) :
            raise ValueError("list_UCD must be a list of UCD medical codes.")

        if (list_CIP13 is not None) and (type(list_CIP13) != list) :
            raise ValueError("list_CIP13 must be a list of CIP-13 medical codes.")

        if (list_ATC is not None) and (type(list_ATC) != list) :
            raise ValueError("list_ATC must be a list of ATC medical codes.")

        if (years is None) or (not isinstance(years, list)): 
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")

        # CCAM
        CCAM_DCIR = self.loc_ccam_dcir(df_ID_PATIENT, years=years, list_CCAM=list_CCAM, print_option=True)
        CCAM_PMSI = self.loc_ccam_pmsi(df_ID_PATIENT, years=years, list_CCAM=list_CCAM, print_option=True, dev=dev) 
        df_dcir = CCAM_DCIR[['BEN_IDT_ANO', 'CAM_PRS_IDE', 'EXE_SOI_DTD']].copy()
        df_dcir = df_dcir.rename(columns={
        'CAM_PRS_IDE': 'CODE',
        'EXE_SOI_DTD': 'DATE'
        })
        df_pmsi = CCAM_PMSI[['BEN_IDT_ANO', 'CDC_ACT', 'EXE_SOI_DTD']].copy()
        df_pmsi = df_pmsi.rename(columns={
        'CDC_ACT': 'CODE',
        'EXE_SOI_DTD': 'DATE'
        })
        CCAM_concat = pd.concat([df_dcir, df_pmsi], ignore_index=True)
        CCAM_concat.sort_values(by=['BEN_IDT_ANO', 'DATE'], inplace=True)
        df_ccam = CCAM_concat.loc[:, ["BEN_IDT_ANO", "DATE", "CODE"]].copy()
        df_ccam["COD_CCAM"] = df_ccam["CODE"]
        df_ccam["COD_ICD10"], df_ccam["COD_CIP"], df_ccam["COD_UCD"], df_ccam["COD_ATC"] = None, None, None, None
        df_ccam = df_ccam[["BEN_IDT_ANO", "DATE", "COD_CCAM", "COD_ICD10", "COD_CIP", "COD_UCD", "COD_ATC"]]

        # ICD10
        ICD10_PMSI = self.loc_icd10_pmsi(df_ID_PATIENT, years=years, list_ICD10=list_ICD10, print_option=True, dev=dev)
        df_icd10 = pd.DataFrame({
        "BEN_IDT_ANO": ICD10_PMSI["BEN_IDT_ANO"],
        "DATE": ICD10_PMSI["EXE_SOI_DTD"],
        "COD_CCAM": None,
        "COD_ICD10": ICD10_PMSI["DGN_PAL"],
        "COD_CIP": None,
        "COD_UCD": None,
        "COD_ATC": None
        })


        # UCD
        UCD_PMSI = self.loc_ucd_pmsi(df_ID_PATIENT, years=years, list_UCD=list_UCD, print_option=True, dev=dev)
        UCD_DCIR = self.loc_ucd_dcir(df_ID_PATIENT, years=years, list_UCD=list_UCD, print_option=True)
        UCD_concat = pd.concat([UCD_DCIR, UCD_PMSI], ignore_index=True)
        UCD_concat.sort_values(by=['BEN_IDT_ANO', 'EXE_SOI_DTD'], inplace=True)

        df_ucd = pd.DataFrame({
        "BEN_IDT_ANO": UCD_concat["BEN_IDT_ANO"],
        "DATE": UCD_concat["EXE_SOI_DTD"],
        "COD_CCAM": None,
        "COD_ICD10": None,
        "COD_CIP": None,
        "COD_UCD": UCD_concat["UCD_UCD_COD"],
        "COD_ATC": None
        })


        # CIP
        CIP_DCIR = self.loc_cip_dcir(df_ID_PATIENT, years=years, list_CIP13=list_CIP13, print_option=True)
        df_cip = pd.DataFrame({
        "BEN_IDT_ANO": CIP_DCIR["BEN_IDT_ANO"],
        "DATE": CIP_DCIR["EXE_SOI_DTD"],
        "COD_CCAM": None,
        "COD_ICD10": None,
        "COD_CIP": CIP_DCIR["PHA_CIP_C13"],
        "COD_UCD": None,
        "COD_ATC": None
        })

        # ATC
        if list_ATC is not None :
            ATC_PMSI = self.loc_atc_pmsi(list_ATC=list_ATC, df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=True, dev=dev)
            ATC_DCIR = self.loc_atc_dcir(list_ATC=list_ATC, df_ID_PATIENT=df_ID_PATIENT, years=years, print_option=True)

            ATC_concat = pd.concat([ATC_DCIR, ATC_PMSI], ignore_index=True)
            ATC_concat.sort_values(by=['BEN_IDT_ANO', 'EXE_SOI_DTD'], inplace=True)

            df_atc = pd.DataFrame({
            "BEN_IDT_ANO": ATC_concat["BEN_IDT_ANO"],
            "DATE": ATC_concat["EXE_SOI_DTD"],
            "COD_CCAM": None,
            "COD_ICD10": None,
            "COD_CIP": None,
            "COD_UCD": None,
            "COD_ATC": ATC_concat["PHA_ATC_CLA"]
            })


        # Fusion
        df_records = pd.concat([df_ccam, df_icd10, df_ucd, df_cip], ignore_index=True)
        df_records = df_records.sort_values(by=["BEN_IDT_ANO", "DATE"]).reset_index(drop=True)

        if export==True :
            df_records.to_pickle(path+'/Bdd.pkl')

        return df_records



        

