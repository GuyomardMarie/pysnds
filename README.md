# pysdns : A Python Package for Exploring the French National Healthcare Database (SNDS)


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

First download the package using either by downloading or cloning the repository
```
git clone https://github.com/GuyomardMarie/pysnds.git
```

Then use one of the following command :

```
pip install -e .
```

or 

```
ou pip install .
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

One specific method, `Get_ID`, retrieves the unique patient identifiers:

- In the DCIR: combination of `BEN_NIR_PSA` and `BEN_RNG_GEM`,
- In the PMSI: `NIR_ANO_17`, linked to `BEN_NIR_PSA` in the DCIR.

This method returns both the updated unique identifier `BEN_IDT_ANO` (ensuring it reflects the latest update date `BEN_DTE_MAJ`) as well as the original identifier pair for completeness.

Once the target population is identified, another method, `Get_AGE`, computes the age at inclusion based on the first occurrence of a relevant medical code.

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

## References
[1] https://zenodo.org/records/15804909?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjQ0OTNlNTkxLWI1YzItNDBhMy04ZjkxLTVkMDc5ODdkMjE2ZCIsImRhdGEiOnt9LCJyYW5kb20iOiIxMGIzMGNhOGE5MDVmNjVhMGU4YjdjYTk1ZmFiMjg3MiJ9._A1LVcH0HRNvKv9MVxV5Q8cgcb1QF1dgNCUBxnm8a-Psx11CnKQ4-XYULeLw1P0HttnPzzFWs4DJa7Wiefs6bA

[2] **The French Early Breast Cancer Cohort (FRESH): a resource for breast cancer research and evaluations of oncology practices based on the French National Healthcare System Database (SNDS)**, E. Dumas, et al., *Cancers*, Vol.14, No.11, p.2671 (2022).

---

## Contribute
Please contact [Marie Guyomard](https://guyomardmarie.github.io/) for reporting any issues or suggest improvements. Collaboration is very welcome, in order to add any functionalities to improve the package and make easier the exploration of the SNDS.





