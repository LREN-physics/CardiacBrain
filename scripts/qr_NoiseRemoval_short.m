function qr_NoiseRemoval_short(Truekallechoes,Truekfitted_realallechoes,Truekfitted_imagallechoes,CoilSensMask,SPMheader,FolderPosition,TE,Maskkspace,WriteNii)
% This function remove the modeled cardiac-induced noise on subset of
% k-space defined using a mask
% 
% qr_NoiseRemoval_short(Truekallechoes,Truekfitted_realallechoes,Truekfitted_imagallechoes,CoilSensMask,SPMheader,FolderPosition,TE,Maskkspace,WriteNii)
%
%
% Input:
%   Truekallechoes            - Original data
%   Truekfitted_realallechoes - Real part of the MR signal that will be
%                               removed bellow the mask
%   Truekfitted_imagallechoes - Imaginary part of the MR signal that will be
%                               removed bellow the mask
%   CoilSensMask              - Mask defining the region where R2* fitting
%                               need to be performed
%   SPMheader                 - Header to plot the .nii
%   FolderPosition            - Output folder of the function
%   TE                        - Echo time for R2* fitting
%   Maskkspace                - Mask of the k-space region that need to be
%   targeted
%   WriteNii                  - Binary to write extra .nii files 
%
% Output:
%   \Threshkspace.nii         - Image of Maskkspace 
%   \DataR2s_NoCard.mat       - Data after targeting the masked region
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

if nargin<9
    WriteNii=1;
end

mkdir(FolderPosition);

% Save Maskkspace in .nii
SPMheader.fname=[FolderPosition,'\Threshkspace.nii'];
qr_spm_write(SPMheader,Maskkspace(:,:,:,1));

disp('Removing modeled cardiac-fluctuation bellow the masked region')
% Remove Truekfitted_realallechoes and Truekfitted_imagallechoes from
% Truekallechoes bellow the mask
Maskkspace=repmat(Maskkspace,1,1,1,size(Truekallechoes,4),size(Truekallechoes,5));
Truek_NoCard_RI=Truekallechoes-(Truekfitted_realallechoes-mean(Truekfitted_realallechoes,4))-1i*(Truekfitted_imagallechoes-mean(Truekfitted_imagallechoes,4));
Truek_NoCard_RI(Maskkspace==0)=Truekallechoes(Maskkspace==0);
True_NoCard_RI_allechoes=repmat(CoilSensMask,1,1,1,size(Truekallechoes,4),size(Truekallechoes,5)).*ifftnjy(Truek_NoCard_RI,[1 2 3]);
clear Threshkspace Truekfitted_realallechoes Truekfitted_imagallechoes

disp('Fitting of the R2s')
% No card RI
R2s_NoCard_RI=zeros(size(Truekallechoes(:,:,:,:,1)));
R2sres_NoCard_RI=zeros(size(Truekallechoes(:,:,:,:,1)));
for cbin=1:size(Truekallechoes,4)
    disp(['No cardiac noise bin ',num2str(cbin),'/',num2str(size(Truekallechoes,4))])
    [R2s_NoCard_RI(:,:,:,cbin),R2sres_NoCard_RI(:,:,:,cbin),~]=qr_R2fitting(squeeze(abs(True_NoCard_RI_allechoes(:,:,:,cbin,:))),TE,CoilSensMask);
end

if WriteNii==1
    SPMheader.fname=[FolderPosition,'\R2s_std_NoCard.nii'];
    qr_spm_write(SPMheader,std(R2s_NoCard_RI,[],4));
    SPMheader.fname=[FolderPosition,'\R2s_mean_NoCard.nii'];
    qr_spm_write(SPMheader,mean(R2s_NoCard_RI,4));
    SPMheader.fname=[FolderPosition,'\R2sres_std_NoCard.nii'];
    qr_spm_write(SPMheader,std(R2sres_NoCard_RI,[],4));
    SPMheader.fname=[FolderPosition,'\R2sres_mean_NoCard.nii'];
    qr_spm_write(SPMheader,mean(R2sres_NoCard_RI,4));
end

disp('Saving the data')
save([FolderPosition,'\DataR2s_NoCard.mat'],'R2s_NoCard_RI','R2sres_NoCard_RI','-v7.3');
clear R2s_NoCard_RI R2sres_NoCard_RI

end