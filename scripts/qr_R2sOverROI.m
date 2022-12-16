function [R2s,R2sres]=qr_R2sOverROI(FolderPosition,AreaOfInterest,Mask)
% Get the R2*/RMSE mean/SD for each brain regions and masks
%
% [R2s,R2sres]=qr_R2sOverROI(FolderPosition,AreaOfInterest,Mask)
%
% Input:
%   FolderPosition  - Folder where the data are located
%   AreaOfInterest  - from \Mask\AreaOfInterest.mat
%   Mask            - Mask used to remove modeled cardiac noise in k-space
%
% Output:
%   R2s             - Mean/SD of R2* over each ROIs for each masks
%   R2sres          - Mean/SD of RMSE over each ROIs for each masks
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

load([FolderPosition,'\Mask_',num2str(1),'\DataR2s_NoCard.mat']);

R2s_Ref=R2s_NoCard_RI;
STDRef=std(R2s_Ref,[],4);
MeanRef=mean(R2s_Ref,4);

for carea=1:length(AreaOfInterest)
    [R2s{carea}.Ref_var(1),R2s{carea}.Ref_var(2)]=qr_AverageOverROI(STDRef,AreaOfInterest{carea}.Mask,AreaOfInterest{carea}.Position);
    [R2s{carea}.Ref_Mean(1),R2s{carea}.Ref_Mean(2)]=qr_AverageOverROI(MeanRef,AreaOfInterest{carea}.Mask,AreaOfInterest{carea}.Position);
end

for cmask=1:length(Mask.Thresh)
    disp(['Mask ',num2str(cmask),'/',num2str(length(Mask.Thresh))])
    load([FolderPosition,'\Mask_',num2str(cmask),'\DataR2s_NoCard.mat']);
    for carea=1:length(AreaOfInterest)
        % R2* std/mean averaged over roi
        [R2s{carea}.R2sSTD_NoCard_RI(cmask,1),R2s{carea}.R2sSTD_NoCard_RI(cmask,2)]=qr_AverageOverROI(std(R2s_NoCard_RI,[],4),AreaOfInterest{carea}.Mask,AreaOfInterest{carea}.Position);
        [R2s{carea}.R2sSTD_NoCard_RI_Mean(cmask,1),R2s{carea}.R2sSTD_NoCard_RI_Mean(cmask,2)]=qr_AverageOverROI(mean(R2s_NoCard_RI,4),AreaOfInterest{carea}.Mask,AreaOfInterest{carea}.Position);
        
        % R2* res std/mean averaged over roi
        [R2sres{carea}.R2sresSTD_NoCard_RI(cmask,1),R2sres{carea}.RatioNoCard_RI(cmask,2)]=qr_AverageOverROI(std(R2sres_NoCard_RI,[],4),AreaOfInterest{carea}.Mask,AreaOfInterest{carea}.Position);
        [R2sres{carea}.R2sresSTD_NoCard_RI_Mean(cmask,1),R2sres{carea}.RatioNoCard_RI_Mean(cmask,2)]=qr_AverageOverROI(mean(R2sres_NoCard_RI,4),AreaOfInterest{carea}.Mask,AreaOfInterest{carea}.Position);
    end
end

end