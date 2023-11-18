# Forked from Nina's FPPE-panel-pipeline

## General workflow

- Read preprocessing: QC of the raw data, trimming, QC of the trimmed data
- Mapping by BWA
- Variant detection and filtering by Mutect2
- CNV calling by CNVkit with different segmentation methods
- Tumor analysis and CNV correction by PureCN:
    - either clustering of CNVs or
    - minimal re-segmentation

## DAG of pipeline in its current configuration
![workflow/DAG/DAG.png](https://github.com/pedricolino/map-and-CNV-calling-benchmark/blob/main/workflow/DAG/DAG.png?raw=true)

