function [Beta,DataFitted]=qr_sinfitharm(Data,HarmonicNumber,Mask)
% This function is used to fit a real signal with sinusoidal functions
%
% [Beta,DataFitted]=sinfitharm(Data,HarmonicNumber,Mask)
%
% Input:
%   Data            - Signal to fit that is 4D, and where the time dimension
%                     is the 4th
%   HarmonicNumber  - Desired number of harmonics that need to be fitted
%   Mask(optional)  - 3D mask that can be used to make the calculation
%                     faster because it will only be done within the mask
%                     
%
% Output:
%   Beta             - Weights of each harmonics
%   DataFitted       - Resulting fitted data
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

if nargin<2
    HarmonicNumber=1;
end
if nargin<3
    Mask=ones(size(Data(:,:,:,1)));
end

if ismatrix(Data)
    if size(Data,1)>size(Data,2)
        Data=permute(Data,[3 4 2 1]);
    else
        Data=permute(Data,[3 4 1 2]);
    end
end

Beta=zeros(size(Data,1),size(Data,2),size(Data,3),1+2*HarmonicNumber);
DataFitted=zeros(size(Data,1),size(Data,2),size(Data,3),size(Data,4));

basefunction=zeros(size(Data,4),2);
basefunction(:,1)=ones(size(Data,4),1);
for ll=1:HarmonicNumber
    basefunction(:,1+2*ll-1)=cos(linspace(0,(2*pi-2*pi/size(Data,4))*ll,size(Data,4)));
    basefunction(:,1+2*ll)=sin(linspace(0,(2*pi-2*pi/size(Data,4))*ll,size(Data,4)));
end

for ii=1:size(Data,1)
    for jj=1:size(Data,2)
        for kk=1:size(Data,3)
            if ~any(squeeze(Data(ii,jj,kk,:))==0)||(Mask(ii,jj,kk)==1)
                Beta(ii,jj,kk,:)=basefunction(~isnan(Data(ii,jj,kk,:)),:)\squeeze(Data(ii,jj,kk,~isnan(Data(ii,jj,kk,:))));
                DataFitted(ii,jj,kk,:)=Beta(ii,jj,kk,1);
                for ll=1:HarmonicNumber
                    DataFitted(ii,jj,kk,:)=squeeze(DataFitted(ii,jj,kk,:))+squeeze(Beta(ii,jj,kk,1+2*ll-1))*cos(linspace(0,(2*pi-2*pi/size(Data,4))*ll,size(Data,4))')+squeeze(Beta(ii,jj,kk,1+2*ll))*sin(linspace(0,(2*pi-2*pi/size(Data,4))*ll,size(Data,4))');
                end
            end
        end
    end
end

end