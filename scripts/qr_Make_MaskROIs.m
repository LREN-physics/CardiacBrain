function qr_Make_MaskROIs(FolderPosition)
% This function is used to create masks of different brain regions
% 
% qr_Make_MaskROIs(FolderPosition)
%
% Input:
%   FolderPosition  - Folder where the data are located
%   Require:
%   - \Labels.mat
%   - \SPMheader.mat
%   - \BetaValues_image.mat
%   - \R2s\R2smap.mat
%   - \STD_data.mat
%
% Output:
%   Create \Mask\AreaOfInterest.mat, with the masks of different 
%   brain areas
%   1: brainstem
%   2: cerebellum
%   3: whole brain
%   4: vessels
%   5: white matter
%   6: grey matter
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

disp('Making the masks')

load([FolderPosition,'\Labels.mat']);
load([FolderPosition,'\SPMheader.mat']);
load([FolderPosition,'\BetaValues_image.mat']);
load([FolderPosition,'\R2s\R2smap.mat']);
load([FolderPosition,'\STD_data.mat']);
TrueR2sMean=mean(TrueR2s,4);

%% Masks

MaskVessels=zeros(size(MaskBrain.GMWM));
% MaskVessels(abs(Beta1(:,:,:,15))>=0.1*max(max(max(abs(Beta1(:,:,:,15))))))=1;
BetaRef=abs(Beta1(:,:,:,15));
% BetaRef=STD_i(:,:,:,15);
BetaRefBrain=BetaRef(MaskBrain.GMWM==1);
MaskVessels(BetaRef>=mean(BetaRefBrain(:)))=1;
MaskVessels(MaskBrain.GMWM==1)=0;

AreaOfInterest{1}.Name='Brainstem';
AreaOfInterest{2}.Name='Cerebellum';
AreaOfInterest{3}.Name='Whole brain';
AreaOfInterest{4}.Name='Vessels';
AreaOfInterest{5}.Name='WM';
AreaOfInterest{6}.Name='GM';

AreaOfInterest{1}.Position=[35,59,60,61,62];%Brainstem
AreaOfInterest{2}.Position=[71,72,73,38,39,40,41];%Cerebelum
AreaOfInterest{3}.Position=[1];%WM+GM
AreaOfInterest{4}.Position=[1];%Vessels
AreaOfInterest{5}.Position=[1];% WM
AreaOfInterest{6}.Position=[1];% GM

AreaOfInterest{1}.Mask=Labels;
AreaOfInterest{2}.Mask=Labels;
AreaOfInterest{3}.Mask=MaskBrain.GMWM;
AreaOfInterest{4}.Mask=MaskVessels;
AreaOfInterest{5}.Mask=MaskBrain.WM;
AreaOfInterest{6}.Mask=MaskBrain.GM;

% Very strigent mask
for carea=1:length(AreaOfInterest)
    % Remove values of R2s that are out of the WM/GM mask
    if carea~=3&&carea~=4
        AreaOfInterest{carea}.Mask(AreaOfInterest{3}.Mask==0)=0;
    end
    % Remove values of R2s that are too high from the masks
    if carea~=4
        AreaOfInterest{carea}.MaskBinary(TrueR2sMean>40)=0;
    end
end

for carea=1:length(AreaOfInterest)
    AreaOfInterest{carea}.MaskBinary=zeros(size(Labels));
    for cpos=1:length(AreaOfInterest{carea}.Position)
        AreaOfInterest{carea}.MaskBinary(AreaOfInterest{carea}.Mask==AreaOfInterest{carea}.Position(cpos))=1;
    end
    SPMheader.fname=[FolderPosition,'\Mask\Mask_',AreaOfInterest{carea}.Name,'.nii'];
    qr_spm_write(SPMheader,AreaOfInterest{carea}.MaskBinary);
end

%%
mkdir([FolderPosition,'\Mask']);
save([FolderPosition,'\Mask\AreaOfInterest.mat'],'AreaOfInterest')

end