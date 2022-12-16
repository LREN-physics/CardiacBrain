function [Mask,Fillingkspace,Threshold]=qr_Mask_kspace(Data,FillingkspaceNeeded)
% This function is used to create a mask covering "FillingkspaceNeeded"
% percents of the largest values of "Data" using an iterative process
% 
% [Mask,Fillingkspace,Threshold]=qr_Mask_kspace(Data,FillingkspaceNeeded)
%
% Input:
%   Data                - Data used to make the 2D mask
%   FillingkspaceNeeded - Requested coverage of the 2D mask
%
% Output:
%   Mask                - Mask covering "FillingkspaceNeeded" percents
%                         of the largest values of "Data"
%	Fillingkspace       - Real percentage of mask coverage
%	Threshold           - Value thresholding between data in/out of mask
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

MaxVal=max(Data(:));
MinVal=min(Data(:));

Fillingkspace(1)=0;
Threshold(1)=(MaxVal+MinVal)./2;

% figure
Mask=zeros(size(Data));
while (abs(Fillingkspace(end)-FillingkspaceNeeded)>=0.01)
    Mask=zeros(size(Data));
    Mask(Data>=Threshold(end))=1;
    Fillingkspace(end+1)=sum(Mask(:))./size(Mask,1)./size(Mask,2)*100;
    
    if Fillingkspace(end)>=FillingkspaceNeeded
        % Filling too low, we need increase the thresh
        MinVal=Threshold(end);
        Threshold(end+1)=(MaxVal+Threshold(end))./2;
    else
        % Filling too high, we need decrease the thresh
        
        MaxVal=Threshold(end);
        Threshold(end+1)=(Threshold(end)+MinVal)./2;
    end
    
    
%     subplot(1,2,1)
%     hold on
%     scatter(length(Fillingkspace),Fillingkspace(end))
%     subplot(1,2,2)
%     hold on
%     scatter(length(Fillingkspace),Threshold(end))
    if length(Fillingkspace)>1000
        break
    end
end


end