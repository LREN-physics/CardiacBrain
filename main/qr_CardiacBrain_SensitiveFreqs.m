function qr_CardiacBrain_SensitiveFreqs(FolderPosition)
% The third main function of the cardiac-induced noise characterization:
% This function is used to find the cardiac noise sensitive k-space region
% It computes maps of R2*/RMSE, R2*/RMSE SD across the cardiac phase of
% the raw data and modeled cardiac noise only, when removing the
% fundamental and first harmonics of the fourier series on a specific
% region of k-space
% 
% qr_CardiacBrain_SensitiveFreqs(FolderPosition)
%
% Input:
%   FolderPosition  - Folder where the data are located
%
% Requires:
%   - \FittedValues_kspace.mat
%   - \SPMheader.mat
%   - \Mask\AreaOfInterest.mat
%   - \Labels.mat
%   - \Truek.mat
%
% Creates:
%   \SensitiveFreqs\CircleMasks.mat  - Save the circular masks
%   \SensitiveFreqs                  - Cardiac noise detrending on raw data
%   \SensitiveFreqs_fitted           - Cardiac noise detrending on modeled
%                                      cardiac-noise only
%   \Results\SensitiveArea.mat
%
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

disp('__________________________________')
disp('Loading data')
load([FolderPosition,'\FittedValues_kspace.mat']);
load([FolderPosition,'\SPMheader.mat']);
load([FolderPosition,'\Mask\AreaOfInterest.mat'])
load([FolderPosition,'\Labels.mat'])
load([FolderPosition,'\Truek.mat']);

%% Mask in k-space creation
disp('__________________________________')
disp('Making masks of k-space regions within which modeled cardiac-induced noise is removed')

mkdir([FolderPosition,'\SensitiveFreqs']);

Nx=size(Truek,2);
Ny=size(Truek,3);

% Create a set of masks covering the center x% of k-space, from 0 to 100 by
% step of 2
Mask.Thresh=linspace(0,100,51);
Mask.Circle=zeros([size(Truek(:,:,:,1)),length(Mask.Thresh)]);
Mask.RefCircle=-squeeze(sqrt(((linspace(1,Nx,Nx)'-(floor(Nx/2)+1)).^2)./Nx+((linspace(1,Ny,Ny)-(floor(Ny/2)+1)).^2)./Ny));
for cmask=1:length(Mask.Thresh)
    [ThreshkspaceCircle,~,~]=qr_Mask_kspace(Mask.RefCircle,Mask.Thresh(cmask));
    Mask.Circle(:,:,:,cmask)=permute(repmat(ThreshkspaceCircle,1,1,size(Mask.Circle,1)),[3 1 2]);
end

save([FolderPosition,'\SensitiveFreqs\CircleMasks.mat'],'Mask')

%% Removing masked noise

% Remove the modeled noise from the raw data, for each masks
disp('__________________________________')
disp('Removing masked noise, on raw data')
mkdir([FolderPosition,'\SensitiveFreqs'])
for cmask=1:length(Mask.Thresh)
    disp(['Mask ',num2str(cmask),'/',num2str(length(Mask.Thresh))])
    qr_NoiseRemoval_short(Truek,Truekfitted_real,Truekfitted_imag,CoilSensMask,SPMheader,[FolderPosition,'\SensitiveFreqs\Mask_',num2str(cmask)],TE,Mask.Circle(:,:,:,cmask));
end

% Remove the modeled noise from the modeled noise, for each masks
disp('__________________________________')
disp('Removing masked noise, on data with only modeled-cardiac noise')
mkdir([FolderPosition,'\SensitiveFreqs_fitted'])
for cmask=1:length(Mask.Thresh)
    disp(['Mask ',num2str(cmask),'/',num2str(length(Mask.Thresh))])
    qr_NoiseRemoval_short(Truekfitted_real+1i*Truekfitted_imag,Truekfitted_real,Truekfitted_imag,CoilSensMask,SPMheader,[FolderPosition,'\SensitiveFreqs_fitted\Mask_',num2str(cmask)],TE,Mask.Circle(:,:,:,cmask));
end

%% Compiling all results and saving

disp('__________________________________')
disp('Computing ROI-wise relative R2*/RMSE SD reduction')

NROIs=length(AreaOfInterest);
NMasks=length(Mask.Thresh);

R2sSTD=zeros(NMasks,NROIs);
R2sSTD_ref=zeros(NROIs);
R2sSTDf=zeros(NMasks,NROIs);
R2sSTD_reff=zeros(NROIs);

R2sresSTD=zeros(NMasks,NROIs);
R2sresSTD_ref=zeros(NROIs);
R2sresSTDf=zeros(NMasks,NROIs);
R2sresSTD_reff=zeros(NROIs);

% Analyis with raw data
[R2s,R2sres]=qr_R2sOverROI([FolderPosition,'\SensitiveFreqs'],AreaOfInterest,Mask);
% Analyis with modeled cardiac-induced noise
[R2sf,R2sresf]=qr_R2sOverROI([FolderPosition,'\SensitiveFreqs_fitted'],AreaOfInterest,Mask);

for carea=1:NROIs
    R2sSTD(:,carea)=squeeze(R2s{carea}.R2sSTD_NoCard_RI(:,1));
    R2sSTD_ref(carea)=squeeze(R2s{carea}.Ref_var(:,1));
    R2sSTDf(:,carea)=squeeze(R2sf{carea}.R2sSTD_NoCard_RI(:,1));
    R2sSTD_reff(carea)=squeeze(R2sf{carea}.Ref_var(:,1));
    
    R2sresSTD(:,carea)=squeeze(R2sres{carea}.R2sresSTD_NoCard_RI(:,1));
    R2sresSTD_ref(carea)=squeeze(R2sres{carea}.R2sresSTD_NoCard_RI(1,1));
    R2sresSTDf(:,carea)=squeeze(R2sresf{carea}.R2sresSTD_NoCard_RI(:,1));
    R2sresSTD_reff(carea)=squeeze(R2sresf{carea}.R2sresSTD_NoCard_RI(1,1));
end

disp('__________________________________')
disp('Saving data...')
save([FolderPosition,'\SensitiveFreqs\R2sSD_allMasks.mat'],'R2sSTD','R2sSTD_ref','R2sresSTD','R2sresSTD_ref','-v7.3');
save([FolderPosition,'\SensitiveFreqs_fitted\R2sSD_allMasks.mat'],'R2sSTDf','R2sSTD_reff','R2sresSTDf','R2sresSTD_reff','-v7.3');


%% Fitting the effect of k-space edges on R2* variability 

disp('__________________________________')
disp('Getting the cardiac-induced noise sensitive k-space region')

mkdir([FolderPosition,'\Results'])
mkdir([FolderPosition,'\Results\Sensitivek'])

% Noise behavior on the edges of k-space goes like sqrt(Npoints)
FitLim=10; % Consider only the 10 last masked regions
DataFitMask=Mask.Thresh(end-FitLim:end);
DataFitMask=DataFitMask-DataFitMask(1);
reg = [DataFitMask'];
W   = (reg'*reg)\reg';
Fsqrt = @(b,x) b.*sqrt(x(end)-x);  % fitted function

FitWeight=zeros(NROIs,1);
NoiseFit=zeros(length(Mask.Thresh),NROIs);
NoiseFitDiff=zeros(length(Mask.Thresh),NROIs);
ElbowPoint=zeros(NROIs,1);

for carea=1:NROIs
    % Data are linearized for the fit
    DataFit=squeeze(R2sSTDf(:,carea));
    DataFit=DataFit(end-FitLim:end);
    DataFit=flip(DataFit').^2;
    beta=W*DataFit';
    
    % Get the non-linearized weights
    FitWeight(carea)=sqrt(beta);
    NoiseFit(:,carea)=Fsqrt(FitWeight(carea),Mask.Thresh)./R2sSTDf(1,carea);
    NoiseFitDiff(:,carea)=(R2sSTDf(:,carea)'-Fsqrt(FitWeight(carea),Mask.Thresh))./R2sSTDf(1,carea);
    ElbowPoint(carea)=qr_find_elbow(Mask.Thresh,NoiseFitDiff(:,carea));
    
    f1=figure;
    subplot(2,2,1)
    hold on
    plot(100-DataFitMask,DataFit)
    plot(100-DataFitMask,beta*DataFitMask)
    title(['Linearized data, ',AreaOfInterest{carea}.Name])
    legend('Original data','Fitted sqrt')
    xlabel('Mask coverage [%]')
    ylabel('(SD(R2*) [s^{-1}])^2')
    subplot(2,2,2)
    hold on
    plot(Mask.Thresh(end-FitLim:end),squeeze(R2sSTDf(end-FitLim:end,carea)))
    plot(Mask.Thresh(end-FitLim:end),Fsqrt(FitWeight(carea),DataFitMask))
    title(['True data, ',AreaOfInterest{carea}.Name])
    legend('Original data','Fitted sqrt')
    xlabel('Mask coverage [%]')
    ylabel('SD(R2*) [s^{-1}]')
    subplot(2,2,3)
    hold on
    plot(Mask.Thresh,R2sSTDf(:,carea)./R2sSTDf(1,carea))
    plot(Mask.Thresh,NoiseFit(:,carea))
    title(['SD(R2*), ',AreaOfInterest{carea}.Name])
    legend('Original data','Fitted sqrt')
    xlabel('Mask coverage [%]')
    ylabel('Relative SD(R2*) [s^{-1}]')
    subplot(2,2,4)
    hold on
    plot(Mask.Thresh,NoiseFitDiff(:,carea))
    plot([ElbowPoint(carea);ElbowPoint(carea)],[min(NoiseFitDiff(:,carea));max(NoiseFitDiff(:,carea))],'k-')
    title(['Original-Fitted SD(R2*), ',AreaOfInterest{carea}.Name])
    legend('Original-Fitted sqrt',['Elbow point at ',num2str(ElbowPoint(carea))])
    xlabel('Mask coverage [%]')
    ylabel('Relative SD(R2*) [s^{-1}]')
    
    saveas(f1,[FolderPosition,'\Results\Sensitivek\Sensitive_regions_',AreaOfInterest{carea}.Name,'.png'])
    close(f1)
end

disp('__________________________________')
disp('Saving data...')
save([FolderPosition,'\SensitiveFreqs\SensitiveArea.mat'],'ElbowPoint');

end