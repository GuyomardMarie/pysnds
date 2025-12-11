# pysnds : A Python Package for Exploring the French National Healthcare Database (SNDS)


The Python package **`pysnds`** simplifies the identification and interpretation of relevant information in the *Système National des Données de Santé*, and therefore the exploitation of this data for research purposes. 

The package offers two key contributions:

1. **Automated Population Identification** – It navigates the SNDS structure to identify specific populations and their characteristics.
2. **Detection of Targeted Medical Events** – It identifies the occurrence of specific medical events within the SNDS for a given population and determines their occurrence dates, including the first appearance.


The package also provides tools to characterize the breast cancer population, which can be easily adapted and customized for other target populations. All the definition of the therapeutical pathways for Breast Cancer are extracted from the FRESH study [1].

---

## Installation

**Requirements:**
- pandas
- numpy
- sqlite
- dateutil
- pkg_resources
- matplotlib

```
git clone https://github.com/GuyomardMarie/pysnds.git
cd pysnds
pip install -e .
```

---

## Description of the package and its functionalities

The `pysnds` package consists of two main classes:

---

### SNDS_Query

This class explores the complex architecture of the SNDS database to identify the target population and retrieve relevant information for the study. It locates medical codes of interest related to:

- Procedures (CCAM-coded) in both the DCIR and PMSI datasets,
- Established diagnoses (ICD-10-coded) in the PMSI,
- Prescribed treatments encoded in CIP-13 format in the DCIR and in UCD format in the PMSI.

All these functionalities take as input unique patient identifiers, a list of targeted medical codes and an optional period of inclusion.

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



---

### SNDS_Treatment

This class inherits from `SNDS_Query`. Given a list of medical codes, it searches the SNDS to determine whether the identified patients have received specific treatments. It also provides the corresponding dates, including the first date of administration.

---

## Use Example

We provide a complete and [reusable example](https://github.com/GuyomardMarie/Test_Package/tree/main/use_example) demonstrating the package functionalities on a simulated dataset [1], illustrated through a case study focused on the **female breast cancer population**. The simulated dataset is provided as a 'pre-processed' dataset. Only treatments related to Breast Cancer are reported.

As a starting point, we include a comprehensive dictionary of medical codes relevant to breast cancer treatments, extracted from the FRESH study [2]:
```
{"Chemotherapy": {
    "ICD10": ["Z511"],
    "CCAM": ["ZZLF900", "EBLA001"],
    "ATC": ["L01AA01", "L01CD02"],
    "UCD": ["9031083", "9031077"]
  },
  "Radiotherapy": {
    "ICD10": ["Z5100", "Z5101"],
    "CCAM": ["ZZNA002", "ZZNL001", "ZZNL002"]
  }
}
```

For each treatment type (surgery, chemotherapy, radiotherapy, etc.), the dictionary contains the corresponding medical codes, structured with keys referring to coding systems (`ATC`, `CCAM`, `ICD10`, `UCD`, `CIP13`). This format is supported by the `SNDS_Query` class, which uses the dictionary to identify the relevant SNDS tables.

The class `SNDS_BC`, built on top of the core classes, characterizes each patient by providing for example:
- Age at enrollment
- Nodal status of the breast cancer
- Type of surgery: none, partial mastectomy, or mastectomy
- Radiotherapy modality: none / neoadjuvant / adjuvant
- Chemotherapy modality: none / neoadjuvant or adjuvant / unitherapy or bitherapy
- Targeted therapy modality: none / neoadjuvant or adjuvant
- Endocrine therapy modality: none / neoadjuvant or adjuvant / detailed regimen
- Breast Cancer Subtype

It also provides tools to make statistical analyzes and visualization tools.

---

## Useful links

- [Health Data Hub Documentation](https://documentation-snds.health-data-hub.fr/snds/tables/)
- [ICD-10](https://icd.who.int/browse10/2019/en)
- [CCAM](https://www.ameli.fr/accueil-de-la-ccam/trouver-un-acte/par-code.php)
- [ATC & CIP13](https://www.vidal.fr/medicaments.html)


---

## References
[1] https://zenodo.org/records/15804909?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjQ0OTNlNTkxLWI1YzItNDBhMy04ZjkxLTVkMDc5ODdkMjE2ZCIsImRhdGEiOnt9LCJyYW5kb20iOiIxMGIzMGNhOGE5MDVmNjVhMGU4YjdjYTk1ZmFiMjg3MiJ9._A1LVcH0HRNvKv9MVxV5Q8cgcb1QF1dgNCUBxnm8a-Psx11CnKQ4-XYULeLw1P0HttnPzzFWs4DJa7Wiefs6bA

[2] **The French Early Breast Cancer Cohort (FRESH): a resource for breast cancer research and evaluations of oncology practices based on the French National Healthcare System Database (SNDS)**, E. Dumas, et al., *Cancers*, Vol.14, No.11, p.2671 (2022).

---

## Contribute
Contributions, issue reports, and feature suggestions are welcome. Please open an issue or a pull request on GitHub, or contact [Marie Guyomard](https://guyomardmarie.github.io/).





