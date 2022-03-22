### Code for:
# Copy Number Variant Risk Scores Associated with Cognition, Psychopathology, and Brain Structure in the Philadelphia Neurodevelopmental Cohort

#### Alexander-Bloch et al., 2022.

This repository contains code to perform CNV quality control, filtering, and annotation as described in this paper. These scripts are performed on output after running CNVision and filtering steps described at: https://github.com/MartineauJeanLouis/MIND-GENESPARALLELCNV

## 1. CNV Processing/Exclusions

CNV exclusion scripts with clean raw CNV data.

| Script | Description |
| --- | --- |
| Exclusion_1 | Takes data from CNVision merge format and applies chip quality control and markers criteria (detailed in SFigs. 1-2).|
| Exclusion_2 | Takes data from Exclusion_1 and applies remaining CNV cleaning criteria (detailed in SFig. 1).|
  
## 2. CNV Annotations

A list of scripts for implementing different CNV annotations, to be run following CNV Exclusions.

| Script | Description |
| --- | --- |
| Annotation_1 | Annotates CNV for pTS and pHI score, returns subject-level data frame for main figures. |
| Annotation_2 | Annotates CNV duplications for main Table 2 & 3. |
| Annotation_3 | Annotates CNV deletions for main Table 2 & 3. |
| Annotation_4 | Annotates CNV for deletion clusters in Figure 1. |
| Annotation_5 | Annotates CNV for duplication clusters in Figure 1. |

## Authors

* Nick Huffnagle
* Aaron Alexander-Bloch
