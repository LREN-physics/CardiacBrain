function [DataAvr,DataSTD]=qr_AverageOverROI(Data,Labels,ROI)
% Average Data over a give region of interest
%
% [R2s,R2sres]=qr_R2sOverROI(FolderPosition,AreaOfInterest,Mask)
%
% Input:
%   Data    - Data that you want to average over the ROI
%   Labels  - Mask with different values for each ROI
%   ROI     - Vector with the different label values of the ROI
%
% Output:
%   DataAvr - Mean value of Data across ROI
%   DataSTD - SD of Data across ROI
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

%% New way

Mask=zeros(size(Labels));
for cROI=1:length(ROI)
    Mask(Labels==ROI(cROI))=1;
end
Mask(isnan(Data))=0;
Mask(isinf(Data))=0;

DataAvr=mean(Data(Mask==1));
DataSTD=std(Data(Mask==1));

end