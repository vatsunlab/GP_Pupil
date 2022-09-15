# GP_Pupil

All code is provided as is.

The parts of the code that extract pupil and motion traces from the data files are hardware specific. The hardware used in our lab (Ripple Neuro) uses the Neuroshare file format. The Matlab Neuroshare API, on which the other code depends, is here: https://github.com/G-Node/nsmatlab. Ensure that the Neuroshare API functions are in the current Matlab path before proceeding with the steps below.


1. Perform a single session pupil analysis first by running pupil_avg_JOVE.m (calls Pupil_analysis_JOVE.m). An example data file (Ring_0001.{ns2, nev, mat}) from a single session is provided. Note that this data file does not replicate the average over sessions and animals plotted in Fig. 3A.

>> pupil_avg_JOVE('GP0000',40,5); %then choose Ring_0001.ns2

2. Combine data across sessions and subjects and run growth curve analysis. To do so, vertically concatenate all the datamat outputs from pupil_avg_JOVE for all sessions, animals,SNRs and attenuations to form the final matrix with [animalID, SNR, dbAtt, Pupil(1-50)timebinvalues]. Use this matrix as input to run the pupil_LME_JOVE code. An example file is provided for the analysis performed in Fig. 3D.

>> load('LME_datamat.mat')
>> pupil_LME_JOVE(LME_datamat)


3. To evaluate reliability of responses, generate another matrix. Put all the session-wise percentage of trials with significant pupil changes into each cell of a cell array where the cells are arranged from lower to higher SNR. Use this cell array as input to run pupil_threshold_estimate_JOVE code. An example file is provided for the analysis performed in Figs. 3B and C.

>> load('pupil_dia_pct_cell_array.mat') % example file for Fig. 3
>> pupil_threshold_estimate_JOVE(pupil_dia_pct_cell_array)
