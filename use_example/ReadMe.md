The Python package *pysnds* simplifies the identification and interpretation of relevant information in the Système National des Données de Santé, and therefore the exploitation of this data for research purposes.

The package offers two key contributions:

- Automated Population Identification – It navigates the SNDS structure to identify specific populations and their characteristics.
Detection of Targeted Medical Events – It identifies the occurrence of specific medical events within the SNDS for a given population and determines their occurrence dates, including the first appearance.
- The package also provides tools to characterize the breast cancer population, which can be easily adapted and customized for other target populations. All the definition of the therapeutical pathways for Breast Cancer are extracted from the FRESH study [1].


**Simulated Data**

The Data used in this tutorial can be found at: https://zenodo.org/records/15804909?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjQ0OTNlNTkxLWI1YzItNDBhMy04ZjkxLTVkMDc5ODdkMjE2ZCIsImRhdGEiOnt9LCJyYW5kb20iOiIxMGIzMGNhOGE5MDVmNjVhMGU4YjdjYTk1ZmFiMjg3MiJ9._A1LVcH0HRNvKv9MVxV5Q8cgcb1QF1dgNCUBxnm8a-Psx11CnKQ4-XYULeLw1P0HttnPzzFWs4DJa7Wiefs6bA

**Table of Contents of the tutorial**

1. Demonstration of the general functionalities of the package

   1.1 Identification of the patients regarding a CCAM code
   
   1.2 Identification of the patients regarding a CIP or UCD code
   
   1.3 Identification of treatment dates
   
   1.4 Identification of the first date treatment
   
   1.5 Search of records of patients
   
   1.6 Identification of the afe at first enrollment

2. Demonstration of the functionalities provided for Breast Cancer analysis
   
     2.1 Get all patients

     2.2 Caracterization of the patients
   
     2.3 Age
   
     2.4 Surgery
   
     2.5 Chemotherapy
   
     2.6 Radiotherapy
   
     2.7 Targeted Therapy
   
     2.8 Endoctrine Therapy
   
     2.9 Nodal Status
   
     2.10 General Function to caracterize the whole population
   
     2.11 Therapeutical Pathways
   
     2.12 Breast Cancer Subtypes
   
     2.13 Statistical Analysis Tools
   
     2.14 Visualization Tools

**Références**

[1] The French Early Breast Cancer Cohort (FRESH): a resource for breast cancer research and evaluations of oncology practices based on the French National Healthcare System Database (SNDS), E. Dumas, et al., Cancers, Vol.14, No.11, p.2671 (2022).
