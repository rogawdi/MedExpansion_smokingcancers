# Medicaid expansion cancer analysis — source-preserving split

This repository is a **mechanical split** of the supplied `smoking_disparity_nonmet_3.R` script. Every source line is preserved in its original order and unchanged within `scripts/`. The original source is retained under `archive/`.

## Scope

The original code includes **combined smoking-related/aerodigestive cancer analyses** as well as head-and-neck, lung, and esophageal sections. It should **not** be described as a verified reproduction of the esophageal cancer JAMA Surgery manuscript: that correspondence has not been established. The named segments are based only on the existing section headings.

## Run order

Read or source `scripts/01_...R` through `scripts/09_...R` **in filename order**, in the same R session, after supplying the authorized source data and appropriate working directory. These files are **not stand-alone scripts**. Some code blocks are exploratory, commented out, or depend on objects or local output paths; the original code has not been executed, debugged, or independently validated.

## Data access and security

The source references SEER-derived CSV data in `inputs/`. These are **not included**. Users must obtain authorized access through the SEER program and follow its use agreements. Do not commit SEER patient-level extracts, credentials, local `.RData`/`.Rhistory` files, or generated patient-level output. The script contains a hard-coded local working-directory path; it is retained intentionally to uphold the no-code-changes constraint and requires local configuration to run.

## Integrity

`SPLIT_MANIFEST.json` provides the exact original line ranges and SHA-256 digest. Concatenating the nine scripts in lexical order reproduces the archived original **byte-for-byte**.

## Release status

**Preservation / organizational copy, not a validated public reproducibility package.** Prior to a public release, identify the authoritative manuscript code, inspect any outputs for restrictions, ensure source data access is compliant, and validate that the scripts reproduce the paper's outputs. No statistical model logic was edited in this split.
