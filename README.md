Medicaid Expansion and Outcomes in Smoking-Related Cancers

Overview

This repository contains R scripts developed to investigate associations between statewide Medicaid expansion under the Affordable Care Act and cancer treatment, survival, and healthcare disparities using data from the Surveillance, Epidemiology, and End Results (SEER) Program.

Analyses include esophageal, lung, and head-and-neck malignancies, as well as combined analyses of smoking-related cancers.

Repository Organization

scripts/ — R scripts organized by analysis section.

archive/ — Original analysis script (smoking_disparity_nonmet_3.R).

SPLIT_MANIFEST.json — Documentation of the original source file and script organization.

The scripts preserve the original analytical code without modifying statistical models, variable definitions, or estimation procedures.

Statistical Analyses

The repository includes code for survival analyses, difference-in-differences models, treatment outcomes, healthcare disparities, and exploratory sensitivity analyses.

Not all analyses contained in the repository were included in the final manuscripts. Individual publications should be consulted for the analyses and results ultimately reported.

Data Availability

Data were obtained from the National Cancer Institute's Surveillance, Epidemiology, and End Results (SEER) Program.

The underlying SEER Research Plus datasets are not distributed in this repository. Researchers wishing to reproduce the analyses must independently obtain access to the appropriate SEER data releases and comply with applicable data-use requirements.

SEER data access: https://seer.cancer.gov/data/access.html

Running the Analyses

Scripts are numbered sequentially and should be executed in filename order within a single R session. They require the corresponding SEER-derived input datasets and installation of the R packages referenced in the code.

Some scripts contain exploratory analyses and references to local file paths. Users may need to configure their local environment and input/output directories before execution.

The source code has been preserved without modification; the reorganized scripts have not been independently executed as a complete pipeline.

Reproducibility and Code Availability

The repository provides the original statistical analysis code for transparency and methodological reference.

The individual scripts were separated from the original analysis file without changing their contents. The accompanying manifest documents the source line ranges and file integrity.

For questions regarding the analyses or specific publications, please contact the corresponding authors of the relevant manuscripts.
