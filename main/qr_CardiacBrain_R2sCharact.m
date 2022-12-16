function qr_CardiacBrain_R2sCharact(FolderPosition)
% The second main function of the cardiac-induced noise characterization:
% This function compute maps of R2*/RMSE, R2*/RMSE SD across the cardiac
% phase of the raw data, raw data when the modeled cardiac noise is
% removed, and modeled cardiac noise only
% 
% qr_CardiacBrain_R2sCharact(FolderPosition)
%
% Input:
%   FolderPosition  - Folder where the data are located
% 
% Requires:
%   - \FittedValues_kspace.mat
%   - \Labels.mat
%   - \Truek.mat
%   - \SPMheader.mat
%
% Creates:
%   \SensitiveFreqs\Reference       - Raw data
%   \SensitiveFreqs\Noise_removed   - Raw data without modeled cardiac
%   noise
%   \SensitiveFreqs\OnlyNoise       - Modeled cardiac-noise only
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

disp('__________________________________')
disp('Loading data')
load([FolderPosition,'\FittedValues_kspace.mat']);
load([FolderPosition,'\Labels.mat']);
load([FolderPosition,'\Truek.mat']);
load([FolderPosition,'\SPMheader.mat']);

disp('__________________________________')
disp('Reference R2*/RMSE SD when no noise is removed - \SensitiveFreqs\Reference')
% Reference R2*/RMSE SD when no noise is removed
qr_NoiseRemoval_short(Truek,Truekfitted_real,Truekfitted_imag,CoilSensMask,SPMheader,[FolderPosition,'\SensitiveFreqs\Reference'],TE,zeros(size(Truek(:,:,:,1,1))));

disp('__________________________________')
disp('R2*/RMSE SD when all fitted cardiac noise is removed - \SensitiveFreqs\Noise_removed')
% R2*/RMSE SD when all fitted cardiac noise is removed
qr_NoiseRemoval_short(Truek,Truekfitted_real,Truekfitted_imag,CoilSensMask,SPMheader,[FolderPosition,'\SensitiveFreqs\Noise_removed'],TE,ones(size(Truek(:,:,:,1,1))));

disp('__________________________________')
disp('Reference R2*/RMSE SD when only fitted cardiac noise is present - \SensitiveFreqs\OnlyNoise')
% Reference R2*/RMSE SD when only fitted cardiac noise is present
qr_NoiseRemoval_short(Truekfitted_real+1i*Truekfitted_imag,zeros(size(Truekfitted_real)),zeros(size(Truekfitted_imag)),CoilSensMask,SPMheader,[FolderPosition,'\SensitiveFreqs\OnlyNoise'],TE,zeros(size(Truek(:,:,:,1,1))));

end