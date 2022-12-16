function qr_overlapImages(Data,Background,Mask,Filename,Scaling1,Scaling2,colormap2,RatioAspect,writetitle)
% Make and save a figure overlapping two images
% 
% qr_overlapImages(Data,Background,Mask,Filename,Scaling1,Scaling2,colormap2,RatioAspect)
%
% Input:
%   Data        - Main image that need to be displayed
%   Background  - Background of the main image for better visualisation,
%                 idealy the magnitude of an anatomical image
%   Mask        - On region where Mask==1, Data is shown, otherwise
%                 Background is shown
%   Filename    - Save the image using this name
%   Scaling1    - Minimum of colormap for Data
%   Scaling2    - Maximum of colormap for Data
%   colormap2   - Color of the colormap
%   RatioAspect - Relative size of the voxels, [1 2 1] for sagital, [1 1 1]
%                 for transverse plan
%
% Output:
%   Make a cool image
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland


Data=squeeze(Data);
Background=squeeze(Background);
Mask=squeeze(Mask);

% Show the change of R2s relative to the mean value
Background(Mask==1)=max(abs(Background(:)));
Data(Mask==0)=nan;

if nargin < 8
    RatioAspect=[1 2 1];
end

% Colomap changed to have the masked region white
% myColorMap = hot(256);
% myColorMap(1,:) = 0.5;

f=figure;
set(gcf, 'Position',  [0, 0, 800, 400])
% set(gcf, 'Position',  [0, 0, 10*size(Data,2), 10*size(Data,1)])
ax1 = axes;
ax2 = axes;
linkaxes([ax1,ax2])
im1 = imagesc(ax1,Background);
set(im1,'AlphaData',0.8*(1-Mask));
daspect([1 2 1])
colormap(ax1,'gray')
set(gcf,'color','w');
im2 = imagesc(Data);
set(im2,'AlphaData',Mask);
caxis(ax2,[Scaling1 Scaling2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
%     axis equal;
daspect(ax1,RatioAspect) % use [1 2 1] n case of non isotropic voxels
daspect(ax2,RatioAspect)
% daspect(ax1,[1 1 1])
% daspect(ax2,[1 1 1])
if nargin<7
    colormap(ax2,hot)
else
    colormap(ax2,colormap2)
end
hbar =colorbar;
set(hbar,'Position',[0.92 0.115 0.03 0.8]);
if exist('writetitle','var')
    set(hbar.Title,{'String','Rotation','Position'},{writetitle,0,[-5 -20]})
end
if nargin>=4
    saveas(f,Filename)
end
end