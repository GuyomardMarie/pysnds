'''
Class for navigating and identifying population in the SNDS
'''

import sqlite3
import pandas as pd
import numpy as np


class SNDS_query() :
    """
    Class for navigating and identifying population in the SNDS.    
    """

    def __init__(self, conn):
        '''
        Parameters
        ----------
        conn : sqlite3.Connection
            Connection to the database.
        '''

        super(SNDS_query, self).__init__()
        self.conn = conn


    def dbGetQuery(self, query):
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

        cursor = self.conn.cursor()
        cursor.execute(query)
        columns = [description[0] for description in cursor.description]
        data = cursor.fetchall()
        df = pd.DataFrame(data, columns=columns)
        cursor.close()
        return df

    

    def Get_ID(self):
        '''
        Method for collecting unique identifiers of the population.

        Returns
        -------
        df : DataFrame
            DataFrame containing the unique identifier of the population in a column named 'BEN_IDT_ANO'.
        '''

        query_ID_patient =  """
            WITH MaxDates AS (
            SELECT BEN_IDT_ANO, MAX(BEN_DTE_MAJ) AS MaxDate
            FROM IR_BEN_R
            GROUP BY BEN_IDT_ANO
            )

            SELECT DISTINCT
            A.BEN_IDT_ANO,
            A.BEN_NIR_PSA, 
            A.BEN_RNG_GEM
            FROM IR_BEN_R A

            INNER JOIN MaxDates B
                    ON A.BEN_IDT_ANO = B.BEN_IDT_ANO AND A.BEN_DTE_MAJ = B.MaxDate
        """

        unique_id = self.dbGetQuery(query_ID_patient)
    
        print('We have ' + str(unique_id.shape[0]) + ' distinct identifiers, ie. breast cancer feamale patients in the database.')

        return unique_id
    

    def Get_AGE(self, df_ID_PATIENT):
        '''
        Method for collecting the population's age at the time of enrollment.

        Returns
        -------
        df_age : DataFrame
            DataFrame containing the age at the time of enrollment (column DATE_DIAG) for each patient (BEN_IDT_ANO).
        '''
            
        liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))

        # DCIR
        query_AGE_DCIR = f"""
            WITH MaxDates AS (
            SELECT BEN_IDT_ANO, MAX(BEN_DTE_MAJ) AS MaxDate
            FROM IR_BEN_R
            GROUP BY BEN_IDT_ANO
            )
                
            SELECT  
            A.BEN_NIR_PSA,      
            A.BEN_RNG_GEM,
            C.BEN_IDT_ANO,  
            A.EXE_SOI_DTD,
            B.BEN_NAI_ANN,
            B.BEN_NAI_MOI
            FROM ER_PRS_F A

            INNER JOIN IR_BEN_R B
                ON B.BEN_NIR_PSA = A.BEN_NIR_PSA AND B.BEN_RNG_GEM = A.BEN_RNG_GEM

            INNER JOIN MaxDates C
                ON B.BEN_IDT_ANO = C.BEN_IDT_ANO AND B.BEN_DTE_MAJ = C.MaxDate

            WHERE C.BEN_IDT_ANO IN ({liste_benidtano_str})
        """

        age_dcir = self.dbGetQuery(query_AGE_DCIR)
        age_dcir['EXE_SOI_DTD'] = pd.to_datetime(age_dcir['EXE_SOI_DTD']) 
        age_dcir['BEN_NAI_MOI'] = pd.to_datetime(age_dcir['BEN_NAI_MOI']) 
        age_dcir = age_dcir.loc[age_dcir.groupby('BEN_IDT_ANO')['EXE_SOI_DTD'].idxmin()].copy()

        def compute_age(row):
            birth_date = pd.to_datetime(row['BEN_NAI_MOI'])
            diag_date = pd.to_datetime(row['EXE_SOI_DTD'])
            age = diag_date.year - birth_date.year - ((diag_date.month, diag_date.day) < (birth_date.month, birth_date.day))
            return age

        age_dcir['AGE_DIAG'] = age_dcir.apply(compute_age, axis=1)


        # PMSI
        query_AGE_PMSI = f"""
            WITH MaxDates AS (
                SELECT BEN_IDT_ANO, MAX(BEN_DTE_MAJ) AS MaxDate
                FROM IR_BEN_R
                GROUP BY BEN_IDT_ANO
            )
            
            SELECT DISTINCT
            C.BEN_IDT_ANO, 
            C.BEN_RNG_GEM, 
            C.BEN_NIR_PSA, 
            B.ETA_NUM, 
            B.RSA_NUM, 
            A.EXE_SOI_AMD, 
            A.EXE_SOI_AMF,
            B.AGE_ANN
            FROM T_MCOaaC A

            INNER JOIN T_MCOaaB B
                ON A.ETA_NUM = B.ETA_NUM AND A.RSA_NUM = B.RSA_NUM

            INNER JOIN IR_BEN_R C
                ON A.NIR_ANO_17 = C.BEN_NIR_PSA

            INNER JOIN MaxDates D
                ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate

            WHERE C.BEN_IDT_ANO IN ({liste_benidtano_str})
            
        """

        age_pmsi = self.dbGetQuery(query_AGE_PMSI)
        age_pmsi['EXE_SOI_AMD'] = pd.to_datetime(age_pmsi['EXE_SOI_AMD'])  
        age_pmsi = age_pmsi.loc[age_pmsi.groupby('BEN_IDT_ANO')['EXE_SOI_AMD'].idxmin()].copy()

        # AGE min
        age_pmsi = age_pmsi.rename(columns={'EXE_SOI_AMD': 'EXE_SOI_DTD', 'AGE_ANN': 'AGE_DIAG'})
        cols = ['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM', 'EXE_SOI_DTD', 'AGE_DIAG']
        df_age = pd.concat([age_dcir[cols], age_pmsi[cols]], ignore_index=True)
        df_age = df_age.loc[df_age.groupby('BEN_IDT_ANO')['AGE_DIAG'].idxmin()].copy()[['BEN_IDT_ANO', 'AGE_DIAG']]

        return df_age


    def loc_ccam_dcir(self, df_ID_PATIENT, list_CCAM, print_option=True):
        '''
        Method for gathering specific CCAM data from the targeted population in the DCIR.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population.
        list_CCAM : list
            List of CCAM codes to be retrieved.
        print_option : bool, optional
            If True, prints the number of unique patients identified with at least one of the specified CCAM codes. Default is True.

        Returns
        -------
        df_ccam_dcir : DataFrame
            DataFrame containing CCAM codes ('CAM_PRS_IDE') from the DCIR, 
            along with their corresponding execution dates ('EXE_SOI_DTD', 'EXE_SOI_DTF'),
            for each patient ('BEN_IDT_ANO').
        '''

        liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
        liste_ccam_str = ', '.join(f"'{valeur}'" for valeur in list_CCAM)

        query_CCAM_DCIR = f"""

            WITH MaxDates AS (
            SELECT BEN_IDT_ANO, MAX(BEN_DTE_MAJ) AS MaxDate
            FROM IR_BEN_R
            GROUP BY BEN_IDT_ANO
            )
                
            SELECT  
            D.BEN_IDT_ANO,      
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

            INNER JOIN MaxDates D
                ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate

            WHERE A.CAM_PRS_IDE IN ({liste_ccam_str}) AND D.BEN_IDT_ANO IN ({liste_benidtano_str})
            """

        df_ccam_dcir = self.dbGetQuery(query_CCAM_DCIR)
        
        if print_option==True :
            print(str(len(np.unique(df_ccam_dcir['BEN_IDT_ANO']))) + ' patient identified using CCAM code in the DCIR.')

        return df_ccam_dcir
    

    def loc_ccam_pmsi(self, df_ID_PATIENT, list_CCAM, print_option=True):
        '''
        Method for gathering specific CCAM data from the targeted population in the PMSI.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population in the PMSI.
        list_CCAM : list
            List of CCAM codes to be retrieved.
        print_option : bool, optional
            If True, prints the number of unique patients identified with at least one of the specified CCAM codes. Default is True.

        Returns
        -------
        df_ccam_pmsi : DataFrame
            DataFrame containing CCAM codes ('CDC_ACT') from the PMSI, 
            along with their execution dates ('EXE_SOI_AMD', 'EXE_SOI_AMF', and 'ENT_DAT_DEL'),
            for each patient ('BEN_IDT_ANO') and specific hospital stays ('ETA_NUM', 'RSA_NUM').
        '''
        
        liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
        liste_ccam_str = ', '.join(f"'{valeur}'" for valeur in list_CCAM)

        query_CCAM_PMSI = f"""
            WITH MaxDates AS (
              SELECT BEN_IDT_ANO, MAX(BEN_DTE_MAJ) AS MaxDate
              FROM IR_BEN_R
              GROUP BY BEN_IDT_ANO
            )
            SELECT DISTINCT
              C.BEN_IDT_ANO, 
              C.BEN_RNG_GEM, 
              C.BEN_NIR_PSA, 
              B.ETA_NUM, 
              B.RSA_NUM, 
              B.CDC_ACT, 
              B.ENT_DAT_DEL, 
              A.EXE_SOI_AMD, 
              A.EXE_SOI_AMF
            FROM T_MCOaaA B

            INNER JOIN T_MCOaaC A 
              ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

            INNER JOIN IR_BEN_R C
              ON A.NIR_ANO_17 = C.BEN_NIR_PSA

            INNER JOIN MaxDates D
              ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate
            
            WHERE B.CDC_ACT IN ({liste_ccam_str}) AND C.BEN_IDT_ANO IN ({liste_benidtano_str})
            """
        
        df_ccam_pmsi = self.dbGetQuery(query_CCAM_PMSI)

        if print_option==True:
            print(str(len(np.unique(df_ccam_pmsi['BEN_IDT_ANO']))) + ' patient identified using CCAM code in the PMSI.')
        
        return df_ccam_pmsi
    
    
    
    def loc_icd10_pmsi(self, df_ID_PATIENT, list_ICD10, print_option=True):
        '''
        Method for gathering specific ICD-10 data from the targeted population in the PMSI.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population in the PMSI.
        list_ICD10 : list
            List of ICD-10 codes to be retrieved.
        print_option : bool, optional
            If True, prints the number of unique patients identified with at least one of the specified ICD-10 codes. Default is True.

        Returns
        -------
        df_icd10_pmsi : DataFrame
            DataFrame containing ICD-10 codes ('DGN_PAL', 'DGN_REL', 'ASS_DGN') from the PMSI, 
            along with their execution dates ('EXE_SOI_AMD', 'EXE_SOI_AMF'), 
            for each patient ('BEN_IDT_ANO') and specific hospital stays ('ETA_NUM', 'RSA_NUM').
        '''
        
        liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
        liste_icd10_str = ', '.join(f"'{valeur}'" for valeur in list_ICD10)

        query_ICD10_PMSI = f"""
            WITH MaxDates AS (
                SELECT BEN_IDT_ANO, MAX(BEN_DTE_MAJ) AS MaxDate
                FROM IR_BEN_R
                GROUP BY BEN_IDT_ANO
            )
            
            SELECT DISTINCT
            C.BEN_IDT_ANO, 
            C.BEN_RNG_GEM, 
            C.BEN_NIR_PSA, 
            B.ETA_NUM, 
            B.RSA_NUM, 
            B.DGN_PAL, 
            B.DGN_REL, 
            NULL AS ASS_DGN,
            A.EXE_SOI_AMD, 
            A.EXE_SOI_AMF
            FROM T_MCOaaB B

            INNER JOIN T_MCOaaC A 
            ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

            INNER JOIN IR_BEN_R C
            ON A.NIR_ANO_17 = C.BEN_NIR_PSA

            INNER JOIN MaxDates D
            ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate

            WHERE B.DGN_PAL IN ({liste_icd10_str}) OR B.DGN_REL IN ({liste_icd10_str}) AND C.BEN_IDT_ANO IN ({liste_benidtano_str})
            
            UNION ALL
            
            SELECT 
            C.BEN_IDT_ANO,
            C.BEN_RNG_GEM, 
            C.BEN_NIR_PSA, 
            E.ETA_NUM, 
            E.RSA_NUM, 
            NULL AS DGN_PAL,
            NULL AS DGN_REL,
            E.ASS_DGN, 
            A.EXE_SOI_AMD, 
            A.EXE_SOI_AMF
            FROM T_MCOaaD E

            INNER JOIN T_MCOaaC A 
            ON E.ETA_NUM = A.ETA_NUM AND E.RSA_NUM = A.RSA_NUM

            INNER JOIN IR_BEN_R C
            ON A.NIR_ANO_17 = C.BEN_NIR_PSA

            INNER JOIN MaxDates D
            ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate

            WHERE E.ASS_DGN IN ({liste_icd10_str}) AND C.BEN_IDT_ANO IN ({liste_benidtano_str})
        """
    
        df_icd10_pmsi = self.dbGetQuery(query_ICD10_PMSI)

        if print_option==True :
            print(str(len(np.unique(df_icd10_pmsi['BEN_IDT_ANO']))) + ' patient identified using ICD10 code in the PMSI.')
        
        return df_icd10_pmsi
    

    def loc_ucd_pmsi(self, df_ID_PATIENT, list_UCD, print_option=True):
        '''
        Method for gathering specific UCD data from the targeted population in the PMSI.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population in the PMSI.
        list_UCD : list
            List of UCD codes to be retrieved.
        print_option : bool, optional
            If True, prints the number of unique patients identified with at least one of the specified UCD codes. Default is True.

        Returns
        -------
        df_ucd_pmsi : DataFrame
            DataFrame containing UCD codes ('UCD_UCD_COD', 'UCD_COD') from the PMSI, 
            along with their execution dates ('EXE_SOI_AMD', 'EXE_SOI_AMF' and 'DELAI'), 
            for each patient ('BEN_IDT_ANO') and specific hospital stays ('ETA_NUM', 'RSA_NUM').
        '''
        
        liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
        liste_ucd_str = ', '.join(f"'{valeur}'" for valeur in list_UCD)

        query_UCD_PMSI = f"""

            WITH MaxDates AS (
                SELECT BEN_IDT_ANO, MAX(BEN_DTE_MAJ) AS MaxDate
                FROM IR_BEN_R
                GROUP BY BEN_IDT_ANO
            )
            SELECT DISTINCT
                C.BEN_IDT_ANO, 
                C.BEN_RNG_GEM, 
                C.BEN_NIR_PSA, 
                B.ETA_NUM, 
                B.RSA_NUM, 
                B.UCD_UCD_COD,
                B.UCD_COD,
                B.DELAI, 
                A.EXE_SOI_AMD, 
                A.EXE_SOI_AMF
            FROM T_MCOaaMED B

            INNER JOIN T_MCOaaC A 
                ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

            INNER JOIN IR_BEN_R C
                ON A.NIR_ANO_17 = C.BEN_NIR_PSA

            INNER JOIN MaxDates D
            ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate

            WHERE B.UCD_UCD_COD IN ({liste_ucd_str}) OR B.UCD_COD IN ({liste_ucd_str}) AND C.BEN_IDT_ANO IN ({liste_benidtano_str})
        """
    
        df_ucd_pmsi = self.dbGetQuery(query_UCD_PMSI)

        if print_option==True :
            print(str(len(np.unique(df_ucd_pmsi['BEN_IDT_ANO']))) + ' patient identified using UCD code in the PMSI.')
        
        return df_ucd_pmsi
    

    def loc_cip_dcir(self, df_ID_PATIENT, list_CIP13, print_option=True):
        '''
        Method for gathering specific UCD data from the targeted population in the DCIR.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population.
        list_CIP13 : list
            List of CIP-13 codes to be retrieved.
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
        
        liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
        liste_cip13_str = ', '.join(f"'{valeur}'" for valeur in list_CIP13)

        query_CIP13_DCIR = f"""
            WITH MaxDates AS (
                SELECT BEN_IDT_ANO, MAX(BEN_DTE_MAJ) AS MaxDate
                FROM IR_BEN_R
                GROUP BY BEN_IDT_ANO
                )
            SELECT DISTINCT
                B.BEN_IDT_ANO,
                A.BEN_NIR_PSA, 
                A.BEN_RNG_GEM, 
                A.DCT_ORD_NUM, 
                A.EXE_SOI_DTD,
                A.EXE_SOI_DTF,
                D.PHA_PRS_C13 AS PHA_CIP_C13,
                E.PHA_ATC_CLA,
                E.PHA_ATC_LIB,
                A.FLX_DIS_DTD, 
                A.FLX_EMT_NUM, 
                A.FLX_EMT_ORD, 
                A.FLX_EMT_TYP, 
                A.FLX_TRT_DTD, 
                A.ORG_CLE_NUM, 
                A.PRS_ORD_NUM, 
                A.REM_TYP_AFF
            FROM ER_PRS_F A

            INNER JOIN IR_BEN_R B
                ON (A.BEN_NIR_PSA = B.BEN_NIR_PSA AND A.BEN_RNG_GEM = B.BEN_RNG_GEM)

            INNER JOIN MaxDates C
                    ON B.BEN_IDT_ANO = C.BEN_IDT_ANO AND B.BEN_DTE_MAJ = C.MaxDate

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

            WHERE                            
                D.PHA_PRS_C13 IN ({liste_cip13_str})
                AND B.BEN_IDT_ANO IN ({liste_benidtano_str})
        """
    
        df_cip_dcir = self.dbGetQuery(query_CIP13_DCIR)

        if print_option==True :
            print(str(len(np.unique(df_cip_dcir['BEN_IDT_ANO']))) + ' patient identified using CIP13 code in the DCIR.')
        
        return df_cip_dcir

