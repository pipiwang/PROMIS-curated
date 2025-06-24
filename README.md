
# Localization of clinically significant prostate cancer on multiparametric MR: an open reproducible analysis on the digitalized PROMIS dataset

Welcome to the repository for our open-source **PROMIS dataset**, along with **preprocessing tools** and code for reproducing **diagnostic accuracy** analysis.

## What's included
- **PROMIS dataset**: An open-source fully digitised dataset curated from the PROMIS study, including aligned radiological and histopathological labels.
- **Preprocessing tools**: Tools to prepare the dataset for subsequent analysis or machine learning tasks. 
- **Diagnostic Accuracy Analysis**: Code to reproduce the main results of the statistical analysis presented in the paper.

## PROMIS Dataset 
### Overview
| Item | Description |
| ---- | ----------- |
| Data Modality | T2-weighted, High-b DWI, ADC |
| Data Format | NifTi |
| Annotations | Lesion countours, prostate gland mask |
| Clinical report | Template biopsy report, radiologist readings |

### Download
You can download the dataset [here](https://exampleurl) by running the following commands:

```
example command line 1
example command lines 2
```
### List of clinical report entry description
Q: do we need a table here to explain the xlsx sheets?
<!-- 
## Preprocessing tools -->

## Generating localised zones on the prostate mask
To generate localised zones of different granularity, run the following script:

  ```bash
  python gen_localised_zones.py
  ```

## Diagnostic Accuracy Analysiss
To reproduce the main analysis results presented in the paper:

1. **Specify Configuration**  
   Define all required variables and directory paths in the `config.py` file. This includes paths to the dataset, output directories, and any relevant parameters.

2. **Run Analysis**  
   Execute the main analysis script:

   ```bash
   python localised_analysis.py
   ```
