function qr_spm_write(SPMheader,data)
% This function is used to create a .nii from the data
% SPM12 is required
%
% qr_spm_write(info,data)
%
% Input:
%   SPMheader   - Header used to create the .nii
%   data        - Data used to create the .nii             
%
% Output:
%   Create a .nii at SPMheader.fname
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

SPMheader.mat(2,1)=-4;
SPMheader.mat(3,2)=2;
SPMheader.mat(1,3)=4;
data=flip(permute(data,[2 1 3]),2);
SPMheader.dim=size(data);
SPMheader.dt(1)=64;
spm_write_vol(SPMheader,data);
end