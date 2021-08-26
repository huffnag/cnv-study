### Working title:
# Omnigenic Impact of Copy Number Variants on Cognition and Psychopathology in the Philadelphia Neurodevelopmental Cohort

#### Authors: TBD

This repository contains code to replicate the main PNC-CNV project methods and results.

## Pre-hoc data structures

These files are currently required prior to running scripts: 

File name | R environment name | Variable names |
| ---  | --- | --- |
| full.cnvision_merge.RData | dtcnv | TBD |
| chip.qc.cag.RData.RData | chip.qc.cag | "ChipID" "Freq" "Burden" "lrrsd" "bafsd" "wavef" "cagid" |
| chipmarkers.csv | chip.markers | "chipids" "nmarkers" |

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

## 3. Analyses/Results

Scripts to take annotated CNV data and generate the main results of the paper.

| Script | Description |
| --- | --- |
| Figure_1 | CNV frequency by pHI/pTS score |
| Figure_2 | CNV effects on cognition and psychopathology |
| Figure_3 | Polygenic risk score results |
| Figure_4 | Brain deviance and CNV burden |
| Table_2 | CNV deletion and duplication effects on overall cognitive accuracy, sorted by lowest AIC |
| Table_3 | CNV intolerance score effects on bi-factor psychosis traits, sorted by lowest AIC |

## 4. Supplement

All supplemental data analyses are accomplished with Supplement.R

## Authors

* Nick Huffnagle
* Aaron Alexander-Bloch
