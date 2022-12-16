function [R2s,R2sres,Beta1]=qr_R2fitting(Data,TE,Mask)
% This function used to fit R2s
% 
% [R2s,R2sres,Beta1]=qr_R2fitting(Data,TE,Mask)
%
% Input:
%   Data            - 4d matrix xyz+echoes
%   TE              - Either the first TE or an array of the TE of all echoes
%   Mask(optional)  - A mask can be used to mask out background and make the
%                     calculation faster
%
% Output:
%   R2s             - R2* from the fit
%   R2sres          - Residuals of the fit std(Data-Fit)
%   Beta1(optional) - Other parameter of the fit, can be used to plot the
%   R2* decay
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

% TE
if size(TE,1)==1 && size(TE,2)==1
    TE=linspace(TE,size(Data,4)*TE,size(Data,4)');
end

if size(TE,1)==1
    TE=TE';
end

% Optional mask
if nargin==2
    Mask=ones(size(Data,1),size(Data,2),size(Data,3));
end

R2s=zeros(size(Data,1),size(Data,2),size(Data,3));
Beta1=zeros(size(Data,1),size(Data,2),size(Data,3));
R2sres=zeros(size(Data,1),size(Data,2),size(Data,3));

% Linearization
Data=log(Data).*repmat(Mask,1,1,1,size(Data,4));
reg = [ones(numel(TE),1) TE];
W   = (reg'*reg)\reg';
Fit = @(b,x) b(1).*exp(b(2).*x);  % fitted function

for xx=1:size(Data,1)
    for yy=1:size(Data,2)
        for zz=1:size(Data,3)
            if Data(xx,yy,zz,1)~=0
                betas=W*squeeze(Data(xx,yy,zz,:));
                R2s(xx,yy,zz)=-betas(2)*1e3;
                Beta1(xx,yy,zz)=exp(betas(1));
                Residuals=Fit([exp(betas(1)),betas(2)],TE)-squeeze(exp(Data(xx,yy,zz,:)));
                R2sres(xx,yy,zz)=std(Residuals(:));
            end
        end
    end
end

% Removing potential issues with the fit
R2s(isnan(R2s))=0;
R2sres(isnan(R2sres))=0;

end