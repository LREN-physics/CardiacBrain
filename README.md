# CardiacNoiseR2s

# Authors
Quentin Raynaud, 2022
Laboratory for Neuroimaging Research
Lausanne University Hospital & University of Lausanne, Lausanne, Switzerland
Copyright (C) 2022 Laboratory for Neuroimaging Research

# Content

The scripts used in CardiacNoiseR2s are the ones used to perform the characterization of cardiac-induced noise.
A subset of data used in the original publication for computations of the results (DOI:). These data contain 5D k-space data (kx/ky/kz/cardiac phase/TE) from one of the participants.

# Main functions

qr_NoiseCharacLowRes_Main is the main function that needs to be run by inputting the folder where the subject data are located.
It requires the 5D k-space data (Truek.mat), the labels (Labels.mat) as well as a header information for SPM (SPMheader).
This function will run the main functions in the following order: qr_CardiacNoiseR2s_Analysis, qr_CardiacNoiseR2s_R2sCharact, qr_CardiacNoiseR2s_SensitiveFreqs and qr_CardiacNoiseR2s_ShowResults.

qr_CardiacNoiseR2s_Analysis:
First main function used to modeled the cardiac-induced noise in image and k-space with fundamental and first harmonic of the Fourier series.
It saves the Fourier series weights (betas) as nifties for visualization using SMP12, and as .mat to be loaded with matlab. It also saves the modeled cardiac-induced noise in k-space.

qr_CardiacNoiseR2s_R2sCharact:
Second main function used to compute maps of R2*/RMSE, R2*/RMSE SD across the cardiac phase of the raw data, raw data when the modeled cardiac noise is removed, and modeled cardiac noise only.

qr_CardiacNoiseR2s_SensitiveFreqs:
Third main used to find the cardiac noise sensitive k-space region. It computes maps of R2*/RMSE, R2*/RMSE SD across the cardiac phase of the raw data and modeled cardiac noise only, when removing the fundamental and first harmonics of the fourier series on a specific region of k-space.

qr_CardiacNoiseR2s_ShowResults:
Final main function that load all the relevant results and plot the results. The different images are saved in \Results

# REQUIREMENTS

The scripts require SPM12 to function properly (avaliable here: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/).

