function qr_CardiacBrain_Analysis(FolderPosition)
% The first main function of the cardiac-induced noise characterization:
% This function is used to modeled the cardiac-induced noise in image and
% k-space
% 
% qr_CardiacBrain_Analysis(FolderPosition)
%
% Input:
%   FolderPosition  - Folder where the data are located
% 
% Requires:
%   - \Labels.mat
%   - \Truek.mat
%   - \SPMheader.mat
%
% Creates:
% \True.mat                 - Save true magnetization data in image space
% \FittedValues_kspace.mat  - Save modeled cardiac-induced noise in k-space
% \BetaValues_kspace.mat    - Save the weight of each harmonics in k-space,
%                             for magnitude and phase
% \BetaValues_kspace_RI.mat - Save the weight of each harmonics in k-space,
%                             for real and imaginary parts
% \R2s\R2smap.mat           - Fitted R2*/RMSE
% \image                    - Noise analysis in image space
% \kspace                   - Noise analysis in k-space
% \xx\SENSE                 - Mean value of the data across cardiac phase
% \xx\mag_1                 - Fundamental frequency of magnitude data
% \xx\mag_2                 - First harmonic of magnitude data
% \xx\phase_1               - Fundamental frequency of phase data
% \xx\phase_2               - First harmonic of phase data
% \xx\real_1                - Fundamental frequency of real part of data
% \xx\real_2                - First harmonic of real part of data
% \xx\imag_1                - Fundamental frequency of imaginary part of data
% \xx\imag_2                - First harmonic of imaginary part of data
% \R2s                      - Fitted R2*/RMSE
% \R2s\Masked               - Fitted R2*/RMSE only showed for brain voxels
%
% Notes:
% The weights of the fourier components are save as complex numbers, in
% which the real part correspond to the weight of the cosinus, and the
% imaginary part the weight of the sinus.
%__________________________________________________________________________
% Copyright (C) 2022 Laboratory for Neuroimaging Research
% Written by Q. Raynaud, 2022.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

disp('__________________________________')
disp(['Starting analysis for data in ',FolderPosition])

disp('Loading data...')
load([FolderPosition,'\Labels.mat']);
load([FolderPosition,'\Truek.mat']);
load([FolderPosition,'\SPMheader.mat']);

TE=2.34*linspace(1,size(Truek,5),size(Truek,5));

True=ifftnjy(Truek,[1,2,3]);
save([FolderPosition,'\True.mat'],'True','TE','CoilSensMask');

%% Fitting data image space

disp('__________________________________')
disp('Cardiac-induced noise characterization in image space')

mkdir([FolderPosition,'\image\SENSE']);
mkdir([FolderPosition,'\image\mag_1']);
mkdir([FolderPosition,'\image\mag_2']);
mkdir([FolderPosition,'\image\phase_1']);
mkdir([FolderPosition,'\image\phase_2']);

Beta0=zeros(size(squeeze(True(:,:,:,1,:))));
Beta1=zeros(size(Beta0));
Beta2=zeros(size(Beta0));

Beta0_phase=zeros(size(Beta0));
Beta1_phase=zeros(size(Beta0));
Beta2_phase=zeros(size(Beta0));
STD_i=zeros(size(Beta0));

for cechoes=1:size(Truek,5)
    % STD
    STD_i(:,:,:,cechoes)=std(abs(True(:,:,:,:,cechoes)),[],4);
    SPMheader.fname=[FolderPosition,'\image\SENSE\STD_i_abs',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,STD_i(:,:,:,cechoes));
end

for cechoes=1:size(True,5)
    
    disp(['Noise characterization in image space, echo ',num2str(cechoes)])
    
    [Beta,~]=qr_sinfitharm(abs(True(:,:,:,:,cechoes)),2);
    [Beta_phase,~]=qr_sinfitharm(unwrap(angle(True(:,:,:,:,cechoes)),[],4),2);
    
    Beta0(:,:,:,cechoes)=Beta(:,:,:,1);
    Beta1(:,:,:,cechoes)=Beta(:,:,:,2)+1i*Beta(:,:,:,3);
    Beta2(:,:,:,cechoes)=Beta(:,:,:,4)+1i*Beta(:,:,:,5);
    
    Beta0_phase(:,:,:,cechoes)=Beta_phase(:,:,:,1);
    Beta1_phase(:,:,:,cechoes)=Beta_phase(:,:,:,2)+1i*Beta_phase(:,:,:,3);
    Beta2_phase(:,:,:,cechoes)=Beta_phase(:,:,:,4)+1i*Beta_phase(:,:,:,5);
    
    % SENSE
    SPMheader.fname=[FolderPosition,'\image\SENSE\Mag',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,Beta0(:,:,:,cechoes));
    SPMheader.fname=[FolderPosition,'\image\SENSE\Phase',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,Beta0_phase(:,:,:,cechoes));
    
    
    % Beta plot magnitude
    SPMheader.fname=[FolderPosition,'\image\mag_1\mag_1_abs_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(Beta1(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\image\mag_2\mag_2_abs_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(Beta2(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\image\mag_1\mag_1_phase_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,angle(Beta1(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\image\mag_2\mag_2_phase_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,angle(Beta2(:,:,:,cechoes)));
    
    temp=angle(Beta1(:,:,:,cechoes));
    temp(temp<0)=temp(temp<0)+pi;
    SPMheader.fname=[FolderPosition,'\image\mag_1\mag_1_phase_0toPi_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,temp);
    temp=angle(Beta2(:,:,:,cechoes));
    temp(temp<0)=temp(temp<0)+pi;
    SPMheader.fname=[FolderPosition,'\image\mag_2\mag_2_phase_0toPi_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,temp);
    
    % Beta plot phase
    SPMheader.fname=[FolderPosition,'\image\phase_1\phase_1_abs_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(Beta1_phase(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\image\phase_2\phase_2_abs_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(Beta2_phase(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\image\phase_1\phase_1_phase_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,angle(Beta1_phase(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\image\phase_2\phase_2_phase_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,angle(Beta2_phase(:,:,:,cechoes)));
    
    temp=angle(Beta1_phase(:,:,:,cechoes));
    temp(temp<0)=temp(temp<0)+pi;
    SPMheader.fname=[FolderPosition,'\image\phase_1\phase_1_phase_0toPi_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,temp);
    temp=angle(Beta2_phase(:,:,:,cechoes));
    temp(temp<0)=temp(temp<0)+pi;
    SPMheader.fname=[FolderPosition,'\image\phase_2\phase_2_phase_0toPi_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,temp);
    
end

%% Fitting data k-space

disp('__________________________________')
disp('Cardiac-induced noise characterization in k-space')

mkdir([FolderPosition,'\kspace\SENSE']);
mkdir([FolderPosition,'\kspace\mag_k1']);
mkdir([FolderPosition,'\kspace\mag_k2']);
mkdir([FolderPosition,'\kspace\phase_k1']);
mkdir([FolderPosition,'\kspace\phase_k2']);
mkdir([FolderPosition,'\kspace\real_k1']);
mkdir([FolderPosition,'\kspace\real_k2']);
mkdir([FolderPosition,'\kspace\imag_k1']);
mkdir([FolderPosition,'\kspace\imag_k2']);

Betak0=zeros(size(squeeze(Truek(:,:,:,1,:))));
Betak1=zeros(size(Betak0));
Betak2=zeros(size(Betak0));

Betak0_phase=zeros(size(Betak0));
Betak1_phase=zeros(size(Betak0));
Betak2_phase=zeros(size(Betak0));

Betak0_real=zeros(size(Betak0));
Betak1_real=zeros(size(Betak0));
Betak2_real=zeros(size(Betak0));

Betak0_imag=zeros(size(Betak0));
Betak1_imag=zeros(size(Betak0));
Betak2_imag=zeros(size(Betak0));

STD_k=zeros(size(Betak0));
Truekfitted_real=zeros(size(Truek));
Truekfitted_imag=zeros(size(Truek));

for cechoes=1:size(Truek,5)
    % STD
    STD_k(:,:,:,cechoes)=std(Truek(:,:,:,:,cechoes),[],4);
    SPMheader.fname=[FolderPosition,'\kspace\SENSE\STD_k_abs',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,STD_k(:,:,:,cechoes));
end

for cechoes=1:size(Truek,5)
    
    disp(['Noise characterization in k-space, echo ',num2str(cechoes)])
    
    [Betak,~]=qr_sinfitharm(abs(Truek(:,:,:,:,cechoes)),2);
    [Betak_phase,~]=qr_sinfitharm(unwrap(angle(Truek(:,:,:,:,cechoes)),[],4),2);
    
    [Betak_real,Truekfitted_real(:,:,:,:,cechoes)]=qr_sinfitharm(real(Truek(:,:,:,:,cechoes)),2);
    [Betak_imag,Truekfitted_imag(:,:,:,:,cechoes)]=qr_sinfitharm(imag(Truek(:,:,:,:,cechoes)),2);
    
    Betak0(:,:,:,cechoes)=Betak(:,:,:,1);
    Betak1(:,:,:,cechoes)=Betak(:,:,:,2)+1i*Betak(:,:,:,3);
    Betak2(:,:,:,cechoes)=Betak(:,:,:,4)+1i*Betak(:,:,:,5);
    
    Betak0_phase(:,:,:,cechoes)=Betak_phase(:,:,:,1);
    Betak1_phase(:,:,:,cechoes)=Betak_phase(:,:,:,2)+1i*Betak_phase(:,:,:,3);
    Betak2_phase(:,:,:,cechoes)=Betak_phase(:,:,:,4)+1i*Betak_phase(:,:,:,5);
    
    Betak0_real(:,:,:,cechoes)=Betak_real(:,:,:,1);
    Betak1_real(:,:,:,cechoes)=Betak_real(:,:,:,2)+1i*Betak_real(:,:,:,3);
    Betak2_real(:,:,:,cechoes)=Betak_real(:,:,:,4)+1i*Betak_real(:,:,:,5);
    
    Betak0_imag(:,:,:,cechoes)=Betak_imag(:,:,:,1);
    Betak1_imag(:,:,:,cechoes)=Betak_imag(:,:,:,2)+1i*Betak_imag(:,:,:,3);
    Betak2_imag(:,:,:,cechoes)=Betak_imag(:,:,:,4)+1i*Betak_imag(:,:,:,5);
    
    % SENSE
    SPMheader.fname=[FolderPosition,'\kspace\SENSE\Mag',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,Betak0(:,:,:,cechoes));
    SPMheader.fname=[FolderPosition,'\kspace\SENSE\Phase',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,Betak0_phase(:,:,:,cechoes));
    SPMheader.fname=[FolderPosition,'\kspace\SENSE\Real',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,Betak0_real(:,:,:,cechoes));
    SPMheader.fname=[FolderPosition,'\kspace\SENSE\Imag',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,Betak0_imag(:,:,:,cechoes));
    
    % Beta plot magnitude
    SPMheader.fname=[FolderPosition,'\kspace\mag_k1\mag_k1_abs_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(Betak1(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\kspace\mag_k2\mag_k2_abs_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(Betak2(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\kspace\mag_k1\mag_k1_phase_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,angle(Betak1(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\kspace\mag_k2\mag_k2_phase_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,angle(Betak2(:,:,:,cechoes)));
    
    temp=angle(Betak1(:,:,:,cechoes));
    temp(temp<0)=temp(temp<0)+pi;
    SPMheader.fname=[FolderPosition,'\kspace\mag_k1\mag_k1_phase_0toPi_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,temp);
    temp=angle(Betak2(:,:,:,cechoes));
    temp(temp<0)=temp(temp<0)+pi;
    SPMheader.fname=[FolderPosition,'\kspace\mag_k2\mag_k2_phase_0toPi_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,temp);
    
    % Beta plot phase
    SPMheader.fname=[FolderPosition,'\kspace\phase_k1\phase_k1_abs_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(Betak1_phase(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\kspace\phase_k2\phase_k2_abs_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(Betak2_phase(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\kspace\phase_k1\phase_k1_phase_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,angle(Betak1_phase(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\kspace\phase_k2\phase_k2_phase_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,angle(Betak2_phase(:,:,:,cechoes)));
    
    temp=angle(Betak1_phase(:,:,:,cechoes));
    temp(temp<0)=temp(temp<0)+pi;
    SPMheader.fname=[FolderPosition,'\kspace\phase_k1\phase_k1_phase_0toPi_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,temp);
    temp=angle(Betak2_phase(:,:,:,cechoes));
    temp(temp<0)=temp(temp<0)+pi;
    SPMheader.fname=[FolderPosition,'\kspace\phase_k2\phase_k2_phase_0toPi_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,temp);
    
    % Beta plot real
    SPMheader.fname=[FolderPosition,'\kspace\real_k1\real_k1_abs_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(Betak1_real(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\kspace\real_k2\real_k2_abs_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(Betak2_real(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\kspace\real_k1\real_k1_phase_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,angle(Betak1_real(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\kspace\real_k2\real_k2_phase_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,angle(Betak2_real(:,:,:,cechoes)));
    
    temp=angle(Betak1_real(:,:,:,cechoes));
    temp(temp<0)=temp(temp<0)+pi;
    SPMheader.fname=[FolderPosition,'\kspace\real_k1\real_k1_phase_0toPi_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,temp);
    temp=angle(Betak2_real(:,:,:,cechoes));
    temp(temp<0)=temp(temp<0)+pi;
    SPMheader.fname=[FolderPosition,'\kspace\real_k2\real_k2_phase_0toPi_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,temp);
    
    % Beta plot imag
    SPMheader.fname=[FolderPosition,'\kspace\imag_k1\imag_k1_abs_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(Betak1_imag(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\kspace\imag_k2\imag_k2_abs_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(Betak2_imag(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\kspace\imag_k1\imag_k1_phase_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,angle(Betak1_imag(:,:,:,cechoes)));
    SPMheader.fname=[FolderPosition,'\kspace\imag_k2\imag_k2_phase_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,angle(Betak2_imag(:,:,:,cechoes)));
    
    temp=angle(Betak1_imag(:,:,:,cechoes));
    temp(temp<0)=temp(temp<0)+pi;
    SPMheader.fname=[FolderPosition,'\kspace\imag_k1\imag_k1_phase_0toPi_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,temp);
    temp=angle(Betak2_imag(:,:,:,cechoes));
    temp(temp<0)=temp(temp<0)+pi;
    SPMheader.fname=[FolderPosition,'\kspace\imag_k2\imag_k2_phase_0toPi_e',num2str(cechoes,'%02d'),'.nii'];
    qr_spm_write(SPMheader,temp);
    
end

%% Saving data

disp('__________________________________')
disp('Saving modeled cardiac noise...')
save([FolderPosition,'\FittedValues_kspace.mat'],'Truekfitted_real','Truekfitted_imag','-v7.3');
save([FolderPosition,'\STD_data.mat'],'STD_i','STD_k','-v7.3');

disp('Saving beta values...')
save([FolderPosition,'\BetaValues_kspace.mat'],'Betak0','Betak0_phase','Betak1','Betak1_phase','Betak2','Betak2_phase','-v7.3');
save([FolderPosition,'\BetaValues_kspace_RI.mat'],'Betak0_real','Betak0_imag','Betak1_real','Betak1_imag','Betak2_real','Betak2_imag','-v7.3');
save([FolderPosition,'\BetaValues_image.mat'],'Beta0','Beta0_phase','Beta1','Beta1_phase','Beta2','Beta2_phase','-v7.3');

%% Fitting R2s decay of each echoes

disp('__________________________________')
disp('R2s fitting of the data...')

mkdir([FolderPosition,'\R2s\']);
TrueR2s=zeros(size(True(:,:,:,:,1)));
TrueR2sres=zeros(size(True(:,:,:,:,1)));
for cbin=1:size(Truek,4)
    disp(['R2s fitting bin ',num2str(cbin),'/',num2str(size(Truek,4))])
    [TrueR2s(:,:,:,cbin),TrueR2sres(:,:,:,cbin),~]=qr_R2fitting(squeeze(abs(True(:,:,:,cbin,:))),TE,CoilSensMask);
    TrueR2s(isnan(TrueR2s))=0;
    TrueR2sres(isnan(TrueR2sres))=0;
    SPMheader.fname=[FolderPosition,'\R2s\R2s_b',num2str(cbin,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(TrueR2s(:,:,:,cbin)));
    SPMheader.fname=[FolderPosition,'\R2s\R2sres_b',num2str(cbin,'%02d'),'.nii'];
    qr_spm_write(SPMheader,abs(TrueR2sres(:,:,:,cbin)));
end

TrueR2sSTD=std(TrueR2s,[],4);
TrueR2sresSTD=std(TrueR2sres,[],4);
SPMheader.fname=[FolderPosition,'\R2s\STD_R2.nii'];
qr_spm_write(SPMheader,TrueR2sSTD);
SPMheader.fname=[FolderPosition,'\R2s\STD_R2res.nii'];
qr_spm_write(SPMheader,TrueR2sresSTD);

mkdir([FolderPosition,'\R2s\Masked']);
for cbin=1:size(TrueR2s,4)
    SPMheader.fname=[FolderPosition,'\R2s\Masked\R2s_b',num2str(cbin,'%02d'),'.nii'];
    qr_spm_write(SPMheader,MaskBrain.GMWM.*abs(TrueR2s(:,:,:,cbin)));
    SPMheader.fname=[FolderPosition,'\R2s\Masked\R2sres_b',num2str(cbin,'%02d'),'.nii'];
    qr_spm_write(SPMheader,MaskBrain.GMWM.*abs(TrueR2sres(:,:,:,cbin)));
end
SPMheader.fname=[FolderPosition,'\R2s\Masked\STD_R2.nii'];
qr_spm_write(SPMheader,MaskBrain.GMWM.*TrueR2sSTD);
SPMheader.fname=[FolderPosition,'\R2s\Masked\STD_R2res.nii'];
qr_spm_write(SPMheader,MaskBrain.GMWM.*TrueR2sresSTD);

disp('Saving R2s maps...')
save([FolderPosition,'\R2s\R2smap.mat'],'TrueR2s','TrueR2sres','TrueR2sSTD','TrueR2sresSTD')

%% Anatomical reference for following masks

mkdir([FolderPosition,'\Mask']);
SPMheader.fname=[FolderPosition,'\Mask\Anatomical.nii'];
qr_spm_write(SPMheader,abs(True(:,:,:,1,1)));
SPMheader.fname=[FolderPosition,'\Mask\Labels.nii'];
qr_spm_write(SPMheader,Labels);
SPMheader.fname=[FolderPosition,'\Mask\MaskGM.nii'];
qr_spm_write(SPMheader,MaskBrain.GM);
SPMheader.fname=[FolderPosition,'\Mask\MaskWM.nii'];
qr_spm_write(SPMheader,MaskBrain.WM);
SPMheader.fname=[FolderPosition,'\Mask\MaskGMWM.nii'];
qr_spm_write(SPMheader,MaskBrain.GMWM);

qr_Make_MaskROIs(FolderPosition)

end