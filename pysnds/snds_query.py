'''
Class for navigating and identifying population in the SNDS
'''

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
        conn : sqlite3.Connection
            Connection to the database.
        '''

        super(SNDS_Query, self).__init__()
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
        Method for collecting unique identifiers of the population in IR_BEN_R.

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
    
        print('We have ' + str(unique_id.shape[0]) + ' distinct identifiers, ie. patients in the database.')

        return unique_id
    

    def Get_AGE(self, df_ID_PATIENT, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], dev=False):
        '''
        Method for collecting the population's age at the time of enrollment.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
        years : list
            List of dates (either years as integers or datetime(yyyy, mm, dd)) defining the period during which to search for CCAM codes. By default 1st of January 2020 and 31 of December 2020.
        dev : bool, optional
            If True, this indicates that the study is performed on a simulated dataset, in which the PMSI does not have any tables referring to specific years. Default is False.

        Returns
        -------
        df_age : DataFrame
            DataFrame containing the age at the time of enrollment (column DATE_DIAG) for each patient (BEN_IDT_ANO).
        '''

        if set(df_ID_PATIENT.columns) != {"BEN_IDT_ANO", "BEN_NIR_PSA", "BEN_RNG_GEM"}:
            raise ValueError(f"df_ID_PATIENT must at least contain the following columns : BEN_IDT_ANO, BEN_NIR_PSA, BEN_RNG_GEM")

        liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))

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
        else:
            end = datetime(years[1], 12, 31)
            year_end = str(years[1])[-2:]

        flxmin = datetime(flxmin_year, 1, 1)
        flxmax_date = (end + relativedelta(months=6)).replace(day=1)
        vecflx = []
        current_date = flxmin
        while current_date <= flxmax_date:
            vecflx.append(current_date.strftime("%Y-%m-%d"))
            current_date += relativedelta(months=1)

        # DCIR
        age_dcir = pd.DataFrame(columns=['BEN_NIR_PSA', 'BEN_RNG_GEM', 'BEN_IDT_ANO', 'EXE_SOI_DTD', 'FLX_DIS_DTD', 'BEN_NAI_ANN', 'BEN_NAI_MOI'])
        for flux in vecflx:
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
                A.FLX_DIS_DTD, 
                B.BEN_NAI_ANN,
                B.BEN_NAI_MOI
                FROM ER_PRS_F A

                INNER JOIN IR_BEN_R B
                    ON B.BEN_NIR_PSA = A.BEN_NIR_PSA AND B.BEN_RNG_GEM = A.BEN_RNG_GEM

                INNER JOIN MaxDates C
                    ON B.BEN_IDT_ANO = C.BEN_IDT_ANO AND B.BEN_DTE_MAJ = C.MaxDate

                WHERE A.FLX_DIS_DTD = '{flux}' AND C.BEN_IDT_ANO IN ({liste_benidtano_str})
            """
            df_flux = self.dbGetQuery(query_AGE_DCIR)
            age_dcir = pd.concat([age_dcir, df_flux], ignore_index=True)

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
        if dev==True:
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
                A.EXE_SOI_DTD, 
                A.EXE_SOI_DTF,
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
            age_pmsi['EXE_SOI_DTD'] = pd.to_datetime(age_pmsi['EXE_SOI_DTD'])  
            age_pmsi = age_pmsi.loc[age_pmsi.groupby('BEN_IDT_ANO')['EXE_SOI_DTD'].idxmin()].copy()

        else :
            age_pmsi = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'ETA_NUM', 'RSA_NUM', 'EXE_SOI_DTD', 'EXE_SOI_DTF', 'AGE_ANN'])
            yy = list(range(int(year_deb), int(year_end)+1))
            for year in yy :
                table_nameB = f"T_MCO{year}B"
                table_nameC = f"T_MCO{year}C"
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
                    A.EXE_SOI_DTD, 
                    A.EXE_SOI_DTF,
                    B.AGE_ANN
                    FROM {table_nameC} A

                    INNER JOIN {table_nameB} B
                        ON A.ETA_NUM = B.ETA_NUM AND A.RSA_NUM = B.RSA_NUM

                    INNER JOIN IR_BEN_R C
                        ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                    INNER JOIN MaxDates D
                        ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate

                    WHERE C.BEN_IDT_ANO IN ({liste_benidtano_str})
                """
                df_flux = self.dbGetQuery(query_AGE_PMSI)
                age_pmsi = pd.concat([age_pmsi, df_flux], ignore_index=True)

            age_pmsi['EXE_SOI_DTD'] = pd.to_datetime(age_pmsi['EXE_SOI_DTD'])  
            age_pmsi = age_pmsi.loc[age_pmsi.groupby('BEN_IDT_ANO')['EXE_SOI_DTD'].idxmin()].copy()

        # AGE min
        age_pmsi = age_pmsi.rename(columns={'EXE_SOI_DTD': 'EXE_SOI_DTD', 'AGE_ANN': 'AGE_DIAG'})
        cols = ['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM', 'EXE_SOI_DTD', 'AGE_DIAG']
        df_age = pd.concat([age_dcir[cols], age_pmsi[cols]], ignore_index=True)
        df_age = df_age.loc[df_age.groupby('BEN_IDT_ANO')['AGE_DIAG'].idxmin()].copy()[['BEN_IDT_ANO', 'AGE_DIAG']]

        return df_age


    def loc_ccam_dcir(self, df_ID_PATIENT=None, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], list_CCAM=None, print_option=True):
        '''
        Method for gathering specific CCAM data from the targeted population in the DCIR.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population. If None, no filter is applied on patients identifiers.
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
            end = years[1]
        else:
            end = datetime(years[1], 12, 31)

        flxmin = datetime(flxmin_year, 1, 1)
        flxmax_date = (end + relativedelta(months=6)).replace(day=1)
        vecflx = []
        current_date = flxmin
        while current_date <= flxmax_date:
            vecflx.append(current_date.strftime("%Y-%m-%d"))
            current_date += relativedelta(months=1)

        df_ccam_dcir = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM', 'CAM_PRS_IDE', 'EXE_SOI_DTD', 'EXE_SOI_DTF', 'FLX_DIS_DTD'])

        for flux in vecflx:
            conditions = []
            if list_CCAM:
                liste_ccam_str = ', '.join(f"'{valeur}'" for valeur in list_CCAM)
                conditions.append(f"A.CAM_PRS_IDE IN ({liste_ccam_str})")
            if df_ID_PATIENT is not None and not df_ID_PATIENT.empty:
                liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
                conditions.append(f"D.BEN_IDT_ANO IN ({liste_benidtano_str})")
            conditions.append(f"B.FLX_DIS_DTD = '{flux}'")
            conditions.append(f"B.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")
            # Construction finale du WHERE
            where_condition = "WHERE " + " AND ".join(conditions)

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
                B.EXE_SOI_DTF,
                B.FLX_DIS_DTD
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

                {where_condition}
                """
            df_flux = self.dbGetQuery(query_CCAM_DCIR)
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
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
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
        conditions.append(f"A.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")
        where_condition = "WHERE " + " AND ".join(conditions)
    
        if dev == True :
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
                    A.EXE_SOI_DTD, 
                    A.EXE_SOI_DTF
                FROM T_MCOaaA B

                INNER JOIN T_MCOaaC A 
                    ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                INNER JOIN IR_BEN_R C
                    ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                INNER JOIN MaxDates D
                    ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate
                
                {where_condition} 
            """
            df_ccam_pmsi = self.dbGetQuery(query_CCAM_PMSI)
            df_ccam_pmsi.drop_duplicates(inplace=True)
            df_ccam_pmsi.reset_index(drop=True, inplace=True)

        else :
            df_ccam_pmsi = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'ETA_NUM', 'RSA_NUM', 'CDC_ACT', 'ENT_DAT_DEL', 'EXE_SOI_DTD', 'EXE_SOI_DTF'])
            yy = list(range(int(year_deb), int(year_end)+1))
            for year in yy :
                table_nameA = f"T_MCO{year}A"
                table_nameC = f"T_MCO{year}C"
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
                        A.EXE_SOI_DTD, 
                        A.EXE_SOI_DTF
                    FROM {table_nameA} B

                    INNER JOIN {table_nameC} A 
                        ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                    INNER JOIN IR_BEN_R C
                        ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                    INNER JOIN MaxDates D
                        ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate
                    
                    {where_condition} 
                """
                df_flux = self.dbGetQuery(query_CCAM_PMSI)
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
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
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
                A.EXE_SOI_DTD, 
                A.EXE_SOI_DTF
                FROM T_MCOaaB B

                INNER JOIN T_MCOaaC A 
                ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                INNER JOIN IR_BEN_R C
                ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                INNER JOIN MaxDates D
                ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate

                {where_condition_1} 
                
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
                A.EXE_SOI_DTD, 
                A.EXE_SOI_DTF
                FROM T_MCOaaD E

                INNER JOIN T_MCOaaC A 
                ON E.ETA_NUM = A.ETA_NUM AND E.RSA_NUM = A.RSA_NUM

                INNER JOIN IR_BEN_R C
                ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                INNER JOIN MaxDates D
                ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate

                {where_condition_2}
            """
            df_icd10_pmsi = self.dbGetQuery(query_ICD10_PMSI)
            df_icd10_pmsi.drop_duplicates(inplace=True)
            df_icd10_pmsi.reset_index(drop=True, inplace=True)


        else :
            df_icd10_pmsi = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'ETA_NUM', 'RSA_NUM', 'CDC_ACT', 'ENT_DAT_DEL', 'EXE_SOI_DTD', 'EXE_SOI_DTF'])
            yy = list(range(int(year_deb), int(year_end)+1))
            for year in yy :
                table_nameB = f"T_MCO{year}B"
                table_nameC = f"T_MCO{year}C"
                table_nameD = f"T_MCO{year}D"
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
                    A.EXE_SOI_DTD, 
                    A.EXE_SOI_DTF
                    FROM {table_nameB} B

                    INNER JOIN {table_nameC} A 
                    ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                    INNER JOIN IR_BEN_R C
                    ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                    INNER JOIN MaxDates D
                    ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate

                    {where_condition_1}
                    
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
                    A.EXE_SOI_DTD, 
                    A.EXE_SOI_DTF
                    FROM {table_nameD} E

                    INNER JOIN {table_nameC} A 
                    ON E.ETA_NUM = A.ETA_NUM AND E.RSA_NUM = A.RSA_NUM

                    INNER JOIN IR_BEN_R C
                    ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                    INNER JOIN MaxDates D
                    ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate

                    {where_condition_2}
                """
                df_flux = self.dbGetQuery(query_ICD10_PMSI)
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
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
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
            conditions.append(f"(B.UCD_UCD_COD IN ({liste_ucd_str}) OR B.UCD_COD IN ({liste_ucd_str}))")
        if df_ID_PATIENT is not None and not df_ID_PATIENT.empty:
            liste_benidtano_str = ', '.join(f"'{val}'" for val in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
            conditions.append(f"C.BEN_IDT_ANO IN ({liste_benidtano_str})")
        conditions.append(f"A.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")
        where_condition = "WHERE " + " AND ".join(conditions)        

        if dev == True :
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
                    A.EXE_SOI_DTD, 
                    A.EXE_SOI_DTF
                FROM T_MCOaaMED B

                INNER JOIN T_MCOaaC A 
                    ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                INNER JOIN IR_BEN_R C
                    ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                INNER JOIN MaxDates D
                ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate

                {where_condition}
            """
            df_ucd_pmsi = self.dbGetQuery(query_UCD_PMSI)
            df_ucd_pmsi.drop_duplicates(inplace=True)
            df_ucd_pmsi.reset_index(drop=True, inplace=True)

        else :
            df_ucd_pmsi = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_RNG_GEM', 'BEN_NIR_PSA', 'ETA_NUM', 'RSA_NUM', 'UCD_UCD_COD', 'UCD_COD', 'DELAI', 'EXE_SOI_DTD', 'EXE_SOI_DTF'])
            yy = list(range(int(year_deb), int(year_end)+1))
            for year in yy :
                table_nameB = f"T_MCO{year}B"
                table_nameC = f"T_MCO{year}C"
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
                        A.EXE_SOI_DTD, 
                        A.EXE_SOI_DTF
                    FROM {table_nameB} B

                    INNER JOIN {table_nameC} A 
                        ON B.ETA_NUM = A.ETA_NUM AND B.RSA_NUM = A.RSA_NUM

                    INNER JOIN IR_BEN_R C
                        ON A.NIR_ANO_17 = C.BEN_NIR_PSA

                    INNER JOIN MaxDates D
                    ON C.BEN_IDT_ANO = D.BEN_IDT_ANO AND C.BEN_DTE_MAJ = D.MaxDate

                    {where_condition}
                """
                df_flux = self.dbGetQuery(query_UCD_PMSI)
                df_ucd_pmsi = pd.concat([df_ucd_pmsi, df_flux], ignore_index=True)

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
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population in the PMSI. If None, no filter is applied on patients identifiers.
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
            end = years[1]
        else:
            end = datetime(years[1], 12, 31)

        flxmin = datetime(flxmin_year, 1, 1)
        flxmax_date = (end + relativedelta(months=6)).replace(day=1)
        vecflx = []
        current_date = flxmin
        while current_date <= flxmax_date:
            vecflx.append(current_date.strftime("%Y-%m-%d"))
            current_date += relativedelta(months=1)

        df_cip_dcir = pd.DataFrame(columns=['BEN_IDT_ANO', 'BEN_NIR_PSA', 'BEN_RNG_GEM', 'DCT_ORD_NUM', 'EXE_SOI_DTD', 'EXE_SOI_DTF', 'PHA_CIP_C13', 'PHA_ATC_CLA', 'PHA_ATC_LIB', 'FLX_DIS_DTD', 'FLX_EMT_NUM', 'FLX_EMT_ORD', 'FLX_EMT_TYP', 'FLX_TRT_DTD', 'ORG_CLE_NUM', 'PRS_ORD_NUM', 'REM_TYP_AFF'])

        for flux in vecflx:
            conditions = []
            if list_CIP13:
                liste_cip13_str = ', '.join(f"'{valeur}'" for valeur in list_CIP13)
                conditions.append(f"D.PHA_PRS_C13 IN ({liste_cip13_str})")
            if df_ID_PATIENT is not None and not df_ID_PATIENT.empty:
                liste_benidtano_str = ', '.join(f"'{valeur}'" for valeur in np.unique(df_ID_PATIENT.BEN_IDT_ANO))
                conditions.append(f"B.BEN_IDT_ANO IN ({liste_benidtano_str})")
            conditions.append(f"A.FLX_DIS_DTD = '{flux}'")
            conditions.append(f"A.EXE_SOI_DTD BETWEEN '{deb.strftime('%Y-%m-%d')}' AND '{end.strftime('%Y-%m-%d')}'")
            where_condition = "WHERE " + " AND ".join(conditions)

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

            {where_condition}
            """
            df_flux = self.dbGetQuery(query_CIP13_DCIR)
            df_cip_dcir = pd.concat([df_cip_dcir, df_flux], ignore_index=True)

        df_cip_dcir.drop_duplicates(inplace=True)
        df_cip_dcir.reset_index(drop=True, inplace=True)

        if print_option==True :
            print(str(len(np.unique(df_cip_dcir['BEN_IDT_ANO']))) + ' patients identified using CIP code in the DCIR.')

        return df_cip_dcir              

    
    
    def Get_records(self, df_ID_PATIENT, years=[datetime(2020, 1, 1), datetime(2020, 12, 31)], list_CCAM=None, list_ICD10=None, list_UCD=None, list_CIP13=None, export=False, path='', dev=False):
        '''
        Method for getting patients records in a DataFrame and optionaly exporting it in pickle format.

        Parameters
        ----------
        df_ID_PATIENT : DataFrame
            DataFrame containing the column 'BEN_IDT_ANO', which holds the unique identifiers of the targeted population.
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
    
        if (years is None) or (not isinstance(years, list)): 
            raise ValueError("`years` must be a list containing the start and end dates, either as integers (years) or as datetime objects (e.g., datetime(yyyy, mm, dd)).")
        
        # CCAM
        CCAM_DCIR = self.loc_ccam_dcir(df_ID_PATIENT, years=years, list_CCAM=list_CCAM, print_option=False)
        CCAM_PMSI = self.loc_ccam_pmsi(df_ID_PATIENT, years=years, list_CCAM=list_CCAM, print_option=False, dev=dev) 
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
        df_ccam["COD_ICD10"], df_ccam["COD_CIP"], df_ccam["COD_UCD"] = None, None, None
        df_ccam = df_ccam[["BEN_IDT_ANO", "DATE", "COD_CCAM", "COD_ICD10", "COD_CIP", "COD_UCD"]]

        # ICD10
        ICD10_PMSI = self.loc_icd10_pmsi(df_ID_PATIENT, years=years, list_ICD10=list_ICD10, print_option=False, dev=dev)
        df_icd10 = pd.DataFrame({
            "BEN_IDT_ANO": ICD10_PMSI["BEN_IDT_ANO"],
            "DATE": ICD10_PMSI["EXE_SOI_DTD"],
            "COD_CCAM": None,
            "COD_ICD10": ICD10_PMSI["DGN_PAL"],
            "COD_CIP": None,
            "COD_UCD": None
        })


        # UCD
        UCD_PMSI = self.loc_ucd_pmsi(df_ID_PATIENT, years=years, list_UCD=list_UCD, print_option=False, dev=dev)
        df_ucd = pd.DataFrame({
            "BEN_IDT_ANO": UCD_PMSI["BEN_IDT_ANO"],
            "DATE": UCD_PMSI["EXE_SOI_DTD"],
            "COD_CCAM": None,
            "COD_ICD10": None,
            "COD_CIP": None,
            "COD_UCD": UCD_PMSI["UCD_UCD_COD"]
        })


        # CIP
        CIP_DCIR = self.loc_cip_dcir(df_ID_PATIENT, years=years, list_CIP13=list_CIP13, print_option=False)
        df_cip = pd.DataFrame({
            "BEN_IDT_ANO": CIP_DCIR["BEN_IDT_ANO"],
            "DATE": CIP_DCIR["EXE_SOI_DTD"],
            "COD_CCAM": None,
            "COD_ICD10": None,
            "COD_CIP": CIP_DCIR["PHA_CIP_C13"],
            "COD_UCD": None
        })


        # Fusion
        df_records = pd.concat([df_ccam, df_icd10, df_ucd, df_cip], ignore_index=True)
        df_records = df_records.sort_values(by=["BEN_IDT_ANO", "DATE"]).reset_index(drop=True)

        if export==True :
            df_records.to_pickle(path+'/Bdd.pkl')

        return df_records



        

