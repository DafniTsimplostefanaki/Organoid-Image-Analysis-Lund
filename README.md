# Image Analysis of 3D Breast Cancer Organoids

This repository contains a Python script developed during my research internship at **Lund University** for the analysis of confocal microscopy images of 3D breast cancer organoids.

The script performs basic image segmentation and quantitative analysis to extract morphological features and spatial protein expression patterns within individual organoids.

---

# Experimental Context

The pipeline was used to analyze and compare organoid development across four distinct experimental conditions:

1) HMW BME (High Molecular Weight Basement Membrane Extract)
2) LMW BME (Low Molecular Weight)
3) HMW Collagen
4) LMW Collagen

The study focused on evaluating how both matrix composition (BME vs. Collagen) and molecular weight (HMW vs. LMW) influence phenotypic shifts and organoid circularity.
This comparison allowed for a systematic assessment of how different biochemical and mechanical properties of the extracellular matrix (ECM) affect tumor cell behavior.

---

# What the Script Does

- Segments organoids using the actin channel and standard thresholding methods.
- Defines central and peripheral regions of each organoid using morphological operations.
- Computes morphological metrics (area, perimeter, circularity).
- Measures mean fluorescence intensities of CK8 and CK14 in the center and periphery.
- Exports results to a CSV file for downstream analysis.

---

# Notes

The same script was used for multiple experimental conditions by modifying only the user configuration block (input path, scale factors, output filename).

Raw microscopy data (`.nd2` files) and experiment-specific metadata are not included in this repository.
