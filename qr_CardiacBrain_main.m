function qr_CardiacBrain_main(FolderPosition)
% This function call all the main functions used for the cardiac-induced
% noise characterization
%
% The data used in this code can be downloaded on Zenodo
% (DOI:10.5281/zenodo.7428605)
% The code necessitates SPM12, that can be downloaded here:
% https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
% 
% qr_CardiacBrain_main(FolderPosition)
%
% Input:
%   FolderPosition  - Folder where the data are located
%                     exemple: 'C:\DATA\CardiacBrain'
%
% Creates:
%   All the data used for noise characterization, check the readme for more
%   information
%
% Requires:
%   - \SPMheader.mat
%   - \Labels.mat
%   - \Truek.mat
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

qr_CardiacBrain_Analysis(FolderPosition)
qr_CardiacBrain_R2sCharact(FolderPosition)
qr_CardiacBrain_SensitiveFreqs(FolderPosition)
qr_CardiacBrain_ShowResults(FolderPosition)

end




