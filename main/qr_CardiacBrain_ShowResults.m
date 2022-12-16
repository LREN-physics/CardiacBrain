function qr_CardiacBrain_ShowResults(FolderPosition)
% The final main function of the cardiac-induced noise characterization:
% This function load all the relevant results and plot the results
% 
% qr_CardiacBrain_R2sCharact(FolderPosition)
%
% Input:
%   FolderPosition  - Folder where the data are located
%
% Creates:
%   \Results\NoiseAmplitude - Maps of the betas in k-space, amplitude of
%                             cardiac-induced noise (Figure 2)
%   \Results\R2smaps        - Maps of the R2*/RMSE mean/SD across the
%                             cardiac phase (Figure 3/4)
%   \Results\Sensitivek     - Relative R2* SD reduction as a function of
%                             the number of detrended k-space frequencies
%                             (Figure 5)
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

disp('__________________________________')
disp('Generating images for data visualization')

%% Load data
disp('__________________________________')
disp('Loading the data...')

mkdir([FolderPosition,'\Results']);

load([FolderPosition,'\Labels.mat']);
load([FolderPosition,'\True.mat']);
load([FolderPosition,'\SPMheader.mat']);
load([FolderPosition,'\FittedValues_kspace.mat']);
load([FolderPosition,'\BetaValues_kspace_RI.mat']);
disp('Done.')

NBins=size(True,4);

%% Noise amplitude in k-space
disp('__________________________________')
disp('Amplitude of cardiac-induced noise in k-space')

Truek_mc_real=Truekfitted_real-repmat(mean(Truekfitted_real,4),[1 1 1 NBins 1]);
Truek_mc_imag=Truekfitted_imag-repmat(mean(Truekfitted_imag,4),[1 1 1 NBins 1]);
NoiseLevel=squeeze(sum(sum(sum(sqrt(Truek_mc_real.^2+Truek_mc_imag.^2),3),2),1));
clear Truek_mc_real Truek_mc_imag

mkdir([FolderPosition,'\Results\NoiseAmplitude']);

f2=figure;
set(gcf, 'Position',  [0, 0, 800, 1000])
ax1 = axes;
imagesc(ax1,squeeze(mean(abs(Betak1_real(:,:,:,1)),1)))
colormap hot
set(gcf,'color','w');
set(gca,'fontsize',25)
caxis([0 4e3])
daspect([1 1 1])
ax1.XTick = [];
ax1.YTick = [];
colorbar
title('\beta_{real} echo 1')
saveas(f2,[FolderPosition,'\Results\NoiseAmplitude\Beta_real_e1.png'])
close(f2)

f2=figure;
set(gcf, 'Position',  [0, 0, 800, 1000])
ax1 = axes;
imagesc(ax1,squeeze(mean(abs(Betak1_real(:,:,:,15)),1)))
colormap hot
set(gcf,'color','w');
set(gca,'fontsize',25)
caxis([0 4e3])
daspect([1 1 1])
ax1.XTick = [];
ax1.YTick = [];
colorbar
title('\beta_{real} echo 15')
saveas(f2,[FolderPosition,'\Results\NoiseAmplitude\Beta_real_e15.png'])
close(f2)

f2=figure;
set(gcf, 'Position',  [0, 0, 800, 1000])
ax1 = axes;
imagesc(ax1,squeeze(mean(abs(Betak1_imag(:,:,:,1)),1)))
colormap hot
set(gcf,'color','w');
set(gca,'fontsize',25)
caxis([0 4e3])
daspect([1 1 1])
ax1.XTick = [];
ax1.YTick = [];
colorbar
title('\beta_{imag} echo 1')
saveas(f2,[FolderPosition,'\Results\NoiseAmplitude\Beta_imag_e1.png'])
close(f2)

f2=figure;
set(gcf, 'Position',  [0, 0, 800, 1000])
ax1 = axes;
imagesc(ax1,squeeze(mean(abs(Betak1_imag(:,:,:,15)),1)))
colormap hot
set(gcf,'color','w');
set(gca,'fontsize',25)
caxis([0 4e3])
daspect([1 1 1])
ax1.XTick = [];
ax1.YTick = [];
colorbar
title('\beta_{imag} echo 15')
saveas(f2,[FolderPosition,'\Results\NoiseAmplitude\Beta_imag_e15.png'])
close(f2)

f2=figure;
set(gcf, 'Position',  [0, 0, 100*length(TE)/2, 100*NBins])
imagesc(TE,linspace(0,2*pi*(1-1/NBins),NBins),1e-8*mean(NoiseLevel,3))
xlabel('Echo time [ms]')
ylabel('Cardiac phase [rad]')
colormap hot
set(gcf,'color','w');
set(gca,'fontsize',25)
set(gca,'YTick',linspace(0,2*pi,5)) 
set(gca,'YTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
hcb=colorbar;
daspect([7 1 1])
caxis([0,6])
% yticks(linspace(0,2*pi*(1-1/NBins),NBins))
% yticklabels(linspace(0,2*pi*(1-1/NBins),NBins))
% daspect([length(TE) NBins 1])
title(["Global amplitude of cardiac"," fluctuation in k-space"])
set(hcb.Title,{'String','Rotation','Position'},{'[a.u.]',0,[25 -35]})
saveas(f2,[FolderPosition,'\Results\NoiseAmplitude\NoiseBinEchoes.png'])
close(f2)

%% R2* variability
disp('__________________________________')
disp('R2*/RMSE variability across cardiac phase')

mkdir([FolderPosition,'\Results\R2smaps']);
load([FolderPosition,'\SensitiveFreqs\Reference\DataR2s_NoCard.mat']);
R2s_Ref=R2s_NoCard_RI;
R2sres_Ref=R2sres_NoCard_RI;
load([FolderPosition,'\SensitiveFreqs\Noise_removed\DataR2s_NoCard.mat']);
R2s_Corr=R2s_NoCard_RI;
R2sres_Corr=R2sres_NoCard_RI;
clear R2s_NoCard_RI R2sres_NoCard_RI

Params.SlicetoShow=20;
Params.SliceCut=70;

% R2s
qr_overlapImages(std(R2s_Ref(Params.SliceCut:end,:,Params.SlicetoShow,:),[],4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2s_std_ref.png'],0,2,hot,[1,2,1],'SD(R_2*) [s^{-1}]')
close
qr_overlapImages(std(R2s_Corr(Params.SliceCut:end,:,Params.SlicetoShow,:),[],4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2s_std_corr.png'],0,2,hot,[1,2,1],'SD(R_2*) [s^{-1}]')
close
qr_overlapImages(std(R2s_Corr(Params.SliceCut:end,:,Params.SlicetoShow,:),[],4)-std(R2s_Ref(Params.SliceCut:end,:,Params.SlicetoShow,:),[],4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2s_std__diff.png'],-1,1,jet,[1,2,1],'SD(R_2*) [s^{-1}]')
close
qr_overlapImages(std(R2s_Corr(Params.SliceCut:end,:,Params.SlicetoShow,:),[],4)./std(R2s_Ref(Params.SliceCut:end,:,Params.SlicetoShow,:),[],4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2s_std__ratio.png'],0,2,hot,[1,2,1],'SD(R_2*) [s^{-1}]')
close

qr_overlapImages(mean(R2s_Ref(Params.SliceCut:end,:,Params.SlicetoShow,:),4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2s_ref.png'],0,40,hot,[1,2,1],'R_2* [s^{-1}]')
close
qr_overlapImages(mean(R2s_Corr(Params.SliceCut:end,:,Params.SlicetoShow,:),4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2s_corr.png'],0,40,hot,[1,2,1],'R_2* [s^{-1}]')
close
qr_overlapImages(mean(R2s_Corr(Params.SliceCut:end,:,Params.SlicetoShow,:),4)-mean(R2s_Ref(Params.SliceCut:end,:,Params.SlicetoShow,:),4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2s__diff.png'],-1,1,jet,[1,2,1],'R_2* [s^{-1}]')
close
qr_overlapImages(mean(R2s_Corr(Params.SliceCut:end,:,Params.SlicetoShow,:),4)./mean(R2s_Ref(Params.SliceCut:end,:,Params.SlicetoShow,:),4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2s__ratio.png'],0,2,jet,[1,2,1],'R_2* [s^{-1}]')
close


% R2s res STD
qr_overlapImages(1e-3*std(R2sres_Ref(Params.SliceCut:end,:,Params.SlicetoShow,:),[],4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2sres_std_ref.png'],0,2,hot,[1,2,1],'SD(RMSE) [a.u.]')
close
qr_overlapImages(1e-3*std(R2sres_Corr(Params.SliceCut:end,:,Params.SlicetoShow,:),[],4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2sres_std_corr.png'],0,2,hot,[1,2,1],'SD(RMSE) [a.u.]')
close
qr_overlapImages(1e-3*std(R2sres_Corr(Params.SliceCut:end,:,Params.SlicetoShow,:),[],4)-1e-3*std(R2sres_Ref(Params.SliceCut:end,:,Params.SlicetoShow,:),[],4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2sres_std__diff.png'],-1,1,jet,[1,2,1],'SD(RMSE) [a.u.]')
close
qr_overlapImages(1e-3*std(R2sres_Corr(Params.SliceCut:end,:,Params.SlicetoShow,:),[],4)./(1e-3*std(R2sres_Ref(Params.SliceCut:end,:,Params.SlicetoShow,:),[],4)),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2sres_std__ratio.png'],0,2,jet,[1,2,1],'SD(RMSE) [a.u.]')
close

qr_overlapImages(1e-3*mean(R2sres_Ref(Params.SliceCut:end,:,Params.SlicetoShow,:),4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2sres_ref.png'],0,30,hot,[1,2,1],'RMSE [a.u.]')
close
qr_overlapImages(1e-3*mean(R2sres_Corr(Params.SliceCut:end,:,Params.SlicetoShow,:),4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2sres_corr.png'],0,30,hot,[1,2,1],'RMSE [a.u.]')
close
qr_overlapImages(1e-3*mean(R2sres_Corr(Params.SliceCut:end,:,Params.SlicetoShow,:),4)-1e-3*mean(R2sres_Ref(Params.SliceCut:end,:,Params.SlicetoShow,:),4),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2sres__diff.png'],-1,1,jet,[1,2,1],'RMSE [a.u.]')
close
qr_overlapImages(1e-3*mean(R2sres_Corr(Params.SliceCut:end,:,Params.SlicetoShow,:),4)./(1e-3*mean(R2sres_Ref(Params.SliceCut:end,:,Params.SlicetoShow,:),4)),mean(abs(True(Params.SliceCut:end,:,Params.SlicetoShow,:,1)),4),MaskBrain.GMWM(Params.SliceCut:end,:,Params.SlicetoShow),[FolderPosition,'\Results\R2smaps\R2sres__ratio.png'],0,2,jet,[1,2,1],'RMSE [a.u.]')
close

%% Sensitive region of k-space
disp('__________________________________')
disp('Sensitive k-space area')

mkdir([FolderPosition,'\Results\Sensitivek'])
load([FolderPosition,'\SensitiveFreqs\CircleMasks.mat'])
load([FolderPosition,'\Mask\AreaOfInterest.mat']);
load([FolderPosition,'\SensitiveFreqs\R2sSD_allMasks.mat']);
load([FolderPosition,'\SensitiveFreqs_fitted\R2sSD_allMasks.mat']);

NROIs=length(AreaOfInterest);
NMasks=length(Mask.Thresh);

ROIsToShow=[1,2,3,4];

f1=figure;
set(gcf,'color','w');
set(gcf, 'Position',  [0, 0, 2000, 500])
for carea=1:length(ROIsToShow)
    subplot(1,length(ROIsToShow),carea)
    hold on
    plot(Mask.Thresh(:),100*(R2sSTD(:,ROIsToShow(carea))./R2sSTD(1,ROIsToShow(carea))),'linewidth',2)
    title(AreaOfInterest{ROIsToShow(carea)}.Name)
    ylim([50,100])
    set(gca,'fontsize',20)
    xlabel('Detrended freq [%]')
    ylabel('Relative R2* SD [s^{-1}]')
end
saveas(f1,[FolderPosition,'\Results\Sensitivek\Sensitive_rawdata.png'])
close

f1=figure;
set(gcf,'color','w');
set(gcf, 'Position',  [0, 0, 2000, 500])
for carea=1:length(ROIsToShow)
    subplot(1,length(ROIsToShow),carea)
    hold on
    plot(Mask.Thresh(:),100*(R2sSTDf(:,ROIsToShow(carea))./R2sSTDf(1,ROIsToShow(carea))),'linewidth',2)
    title(AreaOfInterest{ROIsToShow(carea)}.Name)
    ylim([0,100])
    set(gca,'fontsize',20)
    xlabel('Detrended freq [%]')
    ylabel('Relative R2* SD [s^{-1}]')
end
saveas(f1,[FolderPosition,'\Results\Sensitivek\Sensitive_modeledCard.png'])
close

%%
disp('__________________________________')
disp('All done :)')

end