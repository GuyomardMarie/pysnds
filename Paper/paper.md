---
title: 'pysnds: A Python Package for Exploring the French National Healthcare Database (SNDS)'
tags:
  - Python
  - Système National des Données de Santé
  - Electronic Health Records
authors:
  - name: Marie Guyomard
    corresponding: true
    affiliation: 1
    orcid: 0009-0004-1108-9511
  - name: Anne-Déborah Bouhnik
    affiliation: 1
  - name: Raquel Ureña
    corresponding: true
    affiliation: 1
affiliations:
 - name: SESSTIM, Aix-Marseille University, France
   index: 1
date: 18 November 2025
bibliography: Bib.bib

--- 


## Summary

The French Nationwide Healthcare Database (*Système National des Données de Santé*, SNDS) is a medical claims database covering almost the entire French population and is increasingly used for research purposes, not only for public health purposes but also in Artificial Intelligence. However, due to its complex architecture, the high dimensionality of the data, and the abstract and specific coding of medical events, it remains challenging to navigate, especially for users with limited programming expertise, such as clinicians. 

To tackle this challenge, we propose the Python package `pysnds`, which simplifies the identification and interpretation of relevant information in this complex database, and therefore the exploitation of this data for research purposes. 

The package offers two key contributions:

1. **Automated Population Identification** – It navigates into the SNDS structure to identify specific populations and their characteristics, using specific medical codes.
2. **Detection of Targeted Medical Events** – It identifies the occurrence of specific medical events within the SNDS for a given population and determines their occurrence dates, including the first appearance.

Thus, the `pysnds` package enables a comprehensive understanding of the study population. To showcase its wide range of functionalities, the package provides tools to characterize a population of women being treated for breast cancer. This example can be easily adapted and customized for other target populations.

---

## Statement of need

The French Nationwide Healthcare Database (SNDS), initially created for healthcare reimbursement, contains medical data on approximately 60 million insured patients in France [@tuppin2017value]. It encompasses standardized codes, including diagnoses (ICD-10 classification), medical procedures (Classification Commune des Actes Médicaux, CCAM), and prescribed treatments (Anatomical Therapeutic Chemical classification, ATC). 

The information contained in the SNDS comes from three sources:

- Health insurance expenditure data (e.g., reimbursements for drugs withdrawn from the pharmacy or for medical procedures) provided by the CNAM (*Système National d'Information InterRégime de l'Assurance Maladie*).
- Information on hospital admissions, reasons for admission, medical procedures, and diagnoses provided by ATIH (*Programme de Médicalisation des Systèmes d'Information*).
- Medical causes of death reported by CépiDc (INSERM).

The DCIR is the SNDS database that contains healthcare services provided in towns and cities and submitted for reimbursement by the French national health insurance system. The PMSI MCO records medical activity in the acute care hospital sector, covering public and private medicine, surgery, and obstetrics.

The SNDS has become an essential resource for public health research. It enables large-scale epidemiological studies, analysis of patient care pathways, and prediction of health risks using comprehensive real-world data, thus representing a strategic lever for evaluating and improving health policies.

Several studies have demonstrated the relevance of using SNDS data to investigate the prevalence, incidence, and management of chronic diseases. Notable examples include research on:

- Psoriatic arthritis [@pina2021epidemiologic ; @vegas2022long],
- Dementia [@carcaillon2021prevalence],
- Phenylketonuria [@douillard2023health].

The SNDS also proved particularly useful during the COVID-19 pandemic thanks to near real-time access to data, which eliminated the need for time-consuming and costly post-epidemic surveys [@poucineau2022hospital; @semenzato2021antihypertensive].

Beyond disease analysis, the SNDS offers significant potential for identifying determinants of care trajectories and their economic impact. For instance, studies by [@poulalhon2018use] and [@ossima2023end] on palliative care and end-of-life expenditures highlighted how care pathways significantly influence healthcare costs. These findings underscore the value of medico-administrative data for longitudinal patient monitoring and optimization of care delivery.

Nevertheless, the information within the SNDS is distributed across more than **700 interconnected tables** linked through various join variables.

**Figure 1.** Main SNDS tables used in the provided package. (8*) refers to 8 key variables: `DCT_ORD_NUM`, `FLX_DIS_DTD`, `FLX_EMT_NUM`, `FLX_EMT_ORD`, `FLX_TRT_DTD`, `ORG_CLE_NUM`, `PRS_ORD_NUM`, `REM_TYP_AFF`.

![Main SNDS tables used in the provided package. (8*) refers to 8 key variables: DCT_ORD_NUM, FLX_DIS_DTD, FLX_EMT_NUM, FLX_EMT_ORD, FLX_TRT_DTD, ORG_CLE_NUM, PRS_ORD_NUM, REM_TYP_AFF.](Schema_SNDS.png)


This complex architecture makes it particularly challenging to retrieve relevant clinical information and accurately define target populations, due to the vast number of medical codes involved, which are often abstract and difficult to interpret. To address these challenges, we propose an open-source package designed to facilitate the identification of patient study patients within the SNDS architecture.

---

## Description of the package and its functionalities

The `pysnds` package consists of two main classes: `SNDS_Treatment`, which inherits from the base class `SNDS_Query`. The key attribute of these classes is a connection to the databases. Currently, a **SQLite** connection is supported by the package, adapted for use with the simulated dataset[^1]. In addition, the package supports a **Spark Hive** connection, as it was developed and tested within the Health Data Hub Platform[^2]. The package can also be easily extended to support **Oracle** connections once Python becomes available within the CNAM environment.

[^1] : [Link to the data for reviewers. The dataset will be available in open access when the paper related to its construction will be published.](https://osf.io/n6ctf/?view_only=40eafca1cc7c42d7815abf2d57f7bc14)

[^2] [The Health Data Hub is the structure in charge of granting access to SNDS data.](https://www.health-data-hub.fr/)

---

### SNDS_Query

This class explores the complex architecture of the SNDS database to identify the population of interest and retrieve relevant information for a study. All functionalities of the class take as input a set of unique patient identifiers (a DataFrame containing the columns `BEN_IDT_ANO`, `BEN_NIR_PSA`, and `BEN_RNG_GEM`), a list of targeted medical codes, and an optional inclusion period.

The class includes a function to identify twins within the SNDS. Since it is impossible to determine which medical codes in the PMSI are associated with which individual among twins (see Figure 1 — the variable ``BEN_RNG_GEM is not included in the PMSI), it is recommended not to include twins in the study population.

Several functionalities of the class locate medical codes of interest in different databases (see Figure 1 and Table 1) related to:
- Procedures (CCAM-coded) in both the DCIR and PMSI datasets,
- Established diagnoses (ICD-10-coded) in the PMSI,
- Prescribed treatments encoded in ATC format in both the DCIR and PMSI, and their corresponding CIP-13 (DCIR) and UCD (DCIR and PMSI) codes.


**Table 1.** Functions to locate medical codes with information regarding the databases in the SNDS used to retrieve the information and the outputs. ID* refers to the unique identifiers of a patient `BEN_IDT_ANO`, `BEN_NIR_PSA` and `BEN_RNG_GEM`.


| Function          | Table                              | Columns of output Dataframe                                                                                         | Definition                                                                                                           |
|------------------|------------------------------------|---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------|
| **Identification of CCAM Codes** |                                    |                                                                                                                     |                                                                                                                     |
| loc_ccam_dcir     | ER_CAM_F                            | - ID* <br> - CAM_PRS_IDE <br> - EXE_SOI_DTD <br> - EXE_SOI_DTF                                                     | - Identifiers <br> - Code CCAM <br> - Start Date ('YYYY-MM-DD') <br> - End Date ('YYYY-MM-DD')                     |
| loc_ccam_pmsi     | T_MCOaaA                            | - ID* <br> - CDC_ACT <br> - EXE_SOI_DTD <br> - DELAI <br> - EXE_SOI_DTF                                            | - Identifiers <br> - Code CCAM <br> - Start Date ('YYYY-MM-DD') <br> - Delay from Start Date <br> - End Date ('YYYY-MM-DD') |
| **Identification of ICD-10 Codes** |                           |                                                                                                                     |                                                                                                                     |
| loc_icd10_pmsi    | T_MCOaaB & T_MCOaaD                 | - ID* <br> - DGN_PAL <br> - DGN_REL <br> - ASS_DGN <br> - EXE_SOI_DTD <br> - EXE_SOI_DTF                           | - Identifiers <br> - Principal Diagnosis <br> - Related Diagnosis <br> - Associated Diagnosis <br> - Start Date ('YYYY-MM-DD') <br> - End Date ('YYYY-MM-DD') |
| **Identification of UCD Codes** |                             |                                                                                                                     |                                                                                                                     |
| loc_ucd_dcir      | ER_UCD_F                            | - ID* <br> - UCD_UCD_COD <br> - COD_UCD <br> - PHA_ATC_CLA <br> - PHA_ATC_LIB <br> - PHA_ATC_C03 <br> - PHA_ATC_L03 <br> - EXE_SOI_DTD <br> - EXE_SOI_DTF | - Identifiers <br> - UCD codes CHAR(13) <br> - Last 7 digits of UCD_UCD_COD CHAR(7) <br> - Therapeutic class CHAR(13) <br> - Label of therapeutic class <br> - Second-level ATC code CHAR(3) <br> - Second-level ATC label <br> - Start Date ('YYYY-MM-DD') <br> - End Date ('YYYY-MM-DD') |
| loc_ucd_pmsi      | T_MCOaaMED & T_MCOaaFH              | - ID* <br> - UCD_UCD_COD <br> - COD_UCD <br> - PHA_ATC_CLA <br> - PHA_ATC_LIB <br> - PHA_ATC_C03 <br> - PHA_ATC_L03 <br> - EXE_SOI_DTD <br> - EXE_SOI_DTF | - Identifiers <br> - UCD codes CHAR(13) <br> - Last 7 digits of UCD_UCD_COD CHAR(7) <br> - Therapeutic class CHAR(13) <br> - Label of therapeutic class <br> - Second-level ATC code CHAR(3) <br> - Second-level ATC label <br> - Start Date ('YYYY-MM-DD') <br> - End Date ('YYYY-MM-DD') |
| **Identification of CIP-13 Codes** |                           |                                                                                                                     |                                                                                                                     |
| loc_cip_dcir      | ER_PHA_F                            | - ID* <br> - PHA_PRS_C13 as PHA_CIP_13 <br> - PHA_ATC_CLA <br> - PHA_ATC_LIB <br> - PHA_ATC_C03 <br> - PHA_ATC_L03 <br> - EXE_SOI_DTD <br> - EXE_SOI_DTF | - Identifiers <br> - CIP-13 codes CHAR(13) <br> - Therapeutic class CHAR(13) <br> - Label of therapeutic class <br> - Second-level ATC code CHAR(3) <br> - Second-level ATC label <br> - Start Date ('YYYY-MM-DD') <br> - End Date ('YYYY-MM-DD') |
| **Identification of ATC Codes** |                             |                                                                                                                     |                                                                                                                     |
| loc_atc_dcir      | loc_ucd_dcir <br> loc_cip_dcir      |                                                                                                                     |                                                                                                                     |
| loc_atc_pmsi      | loc_ucd_pmsi                        |                                                                                                                     |                                                                                                                     |



Finally, the package provides a function to retrieve medical records containing specific medical codes (ICD-10, CCAM, ATC, UCD, and CIP-13) over a defined period for the target population and exports the results in Pickle format.

---

### SNDS_Treatment

This class inherits from `SNDS_Query`. Given a list of medical codes (CCAM, ICD-10, ATC, UCD and CIP-13) and optionaly a period of inclusion, it searches the SNDS to determine whether the identified patients (a DataFrame containing the columns `BEN_IDT_ANO`, `BEN_NIR_PSA`, and `BEN_RNG_GEM`) have received specific treatments. It also provides the corresponding dates, including the first date of administration.

---

## Use example on Female Breast Cancer Population

We provide a complete and reusable example demonstrating the package functionalities, illustrated through a case study focused on the **female breast cancer population**.

We include a comprehensive dictionary (json format) of medical codes related to breast cancer treatments, extracted from the FRESH study [@dumas2022french] as illustrated in Figure 2 :

**Figure 2.** Illustration of the furnished dictionary.

![Illustration of the furnished dictionary.](BC_medical_codes.png)

For each treatment type (surgery, chemotherapy, radiotherapy, etc.), the dictionary contains the corresponding medical codes, structured with keys referring to coding systems (`CCAM`, `ICD10`, `ATC`, `UCD`, `CIP13`). This format is supported by the `SNDS_Query` class, which uses the dictionary to identify the relevant SNDS tables.

The class `SNDS_BC`, built on top of the core classes, analyzes the target population. It includes:

- **Date_Diag** - Determines the date of diagnosis for each patient.
- **Get_Age** - Determines the age of patients during a certain period.
- **Age_Diagnosis** - Determines the age at the date of diagnosis for each patient.
- **treatment_setting** – Determines whether a treatment is neoadjuvant (before surgery) or adjuvant (after surgery).
- **Chemotherapy_Regimen** – Identifies the chemotherapy regimen (unitherapy or bitherapy).
- **EndoctrineTherapy_Treatment** – Identifies the endocrine therapy regimen. Seven regimens are defined:
  - Two unitherapies (Aromatase Inhibitor or Tamoxifen),
  - Four bitherapies (e.g., Tamoxifen followed by Aromatase Inhibitor),
  - One undefined ("Unknown") when the regimen cannot be classified.

**Figure 3.** Example of Characterization of the simumlated Breast Cancer Population using pysnds.
![Example of Characterization of the Breast Cancer Population using pysnds.](Pop_Desc.png)

Based on these functionalities, the `BC_POP_Stat` function characterizes each patient by providing :

- Age at enrollment
- Nodal status of the breast cancer
- Type of surgery: none, partial mastectomy, or mastectomy
- Radiotherapy modality: none / neoadjuvant / adjuvant
- Chemotherapy modality: none / neoadjuvant or adjuvant / unitherapy or bitherapy
- Targeted therapy modality: none / neoadjuvant or adjuvant
- Endocrine therapy modality: none / neoadjuvant or adjuvant / detailed regimen


Finaly, some additional functions contained in this class detects the breast cancer subtype (**BC_subtype**), the different therapeutic sequence, as 11 pathways are defined (**therapeutic_pathway**). We also provide statistical and visualization tools (**statistical_analyses**  and **vizualisation_pop**).

The Figure 3 is produced by these functionalities and illustrates the different caracterizations provided by this class on the simulated dataset[^3].

[^3] : [Link to the data for reviewers. The dataset will be available in open access when the paper related to its construction will be published.](https://osf.io/n6ctf/?view_only=40eafca1cc7c42d7815abf2d57f7bc14)

Although tailored to breast cancer, these tools can be adapted to other pathologies by users, by providing a dictionnary of medical codes related to the specific pathology. A complete example of use of the whole functionalities is provided with the package.

---

## Conclusion and Perspectives

The provided package allows SNDS users to efficiently access relevant medical information within this complex database architecture. Developed in `Python` for compatibility with the Health Data Hub platform, it can also be adapted to `R` to extend its use to the CNAM platform. While the illustrative example focuses on breast cancer, the approach can be readily applied to other pathologies, provided that the corresponding treatment dictionaries are available. Contributions and collaborations are highly encouraged to further enrich the package, enhance its functionalities, and facilitate the exploration of the SNDS.

---

## Acknowledgements

We thank Thomas Guyet for useful discussions about the architecture of the SNDS.

This work was supported by a research grant from the Research French Agency ANR - PEPR Digital Health under the number ANR-22-PESN-0005.The project leading to this publication has received funding from the Excellence Initiative of Aix-Marseille Université - AMidex, a French “Investissements d’Avenir programme” AMX-21-IET-017.

---

## References

