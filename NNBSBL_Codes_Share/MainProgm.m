clear all
close all
clc
%-----------Load the essential files for source localization, which are stored in the "ParamMat" directory----------%
%Gain_constrained --- Lead field matrix, dim: (number of channels)¡Á(number of voxels)
%Faces --- A matrix whose each row consists of the indexes of 3 vertices forming a triangle
%GridLoc --- A matrix consisting of the 3D coordinates of the voxels, dim: (number of voxels)¡Á3
%Region --- A cell consisting of the abbreviative names of the ROIs
%Vertices1 --- A cell consisting of the indexes of vertices (voxels) of the ROIs
ParamMat_name_data = dir([pwd,'\ParamMat\*.mat']);
ParamMat_dir_data = [pwd,'\ParamMat\'];
for kk = 1:size(ParamMat_name_data,1)
    ParamMatName = ParamMat_name_data(kk,1).name;
    filename = [ParamMat_dir_data ParamMatName];
    load(filename);
end
Gain_constrained([19,40],:)=[];%These two channels are not EEG recording channels, hence being excluded.
%-------------------------------------------------------------------------------------------------------------------------%
%Generate a cell-style "DeltaL", whose elements are Laplace Operator matrices for ROIs, and the columns of each Laplace Operator matrix are arranged as same as those in Vertices1.
%Generate a cell-style "C", where C{k}=pinv(DeltaL{k}'*DeltaL{k}).
[DeltaL,C] = LaplaceOperMatGen(Faces,GridLoc,Vertices1);
rangecore = brainrangecluster(GridLoc,Vertices1);%Each ROI's central location
%Form the basis matrices "A_estimate" for the proposed algorithm
Am_estimate = [];
regionweight = zeros(length(DeltaL),1);
for i = 1:length(DeltaL)
    Vertices_temp = Vertices1{i};
    Atemp = Gain_constrained(:,Vertices_temp)*C{i}*Gain_constrained(:,Vertices_temp)';
    regionweight(i) = norm(Atemp,'fro');
    Atemp = Atemp/regionweight(i);
    A(i,:,:) = Atemp;
end
M = size(Gain_constrained,1);
for m = 1:M
    A_m_temp = [];
    for i = 1:length(DeltaL)
        A_temp = A(i,:,m);
        A_m_temp(:,i) = A_temp;
    end
    Am_estimate(:,:,m) = A_m_temp;
end
clear A A_m_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Construct the matrices "U" and "H" used in compared algorithm BSBL-2S %%%%%%%%%%%%%%%
[U,s2,H,Hg] = BSBL2SMatBuild(Gain_constrained,Vertices1,DeltaL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters for generating 
SampleRate = 1000;%Sampling rate
SNR = 15;%SNR
dattimlen = 0.080;%data duration (s)
Ns = floor(SampleRate*dattimlen);%Number of samples
varsigrange = [1,1];%Random fluctuation range of source signal variances
K = length(Region);%Number of all ROIs
% source_num = 4;%Number of activated ROIs
% ind_temp = randperm(K,source_num);
% SrcExistRegion = sort(ind_temp,'ascend');
SrcExistRegion = [75,82,85,128];%The indexes of activated ROIs
source_num = length(SrcExistRegion);
gamma = zeros(K,1);
%Ganerate the noisy multi-channel EEG data "Y", dim: Nch¡ÁNs, where Nch is the number of channels.
%gamma is a vector contains the variance paramters for the ROIs
[Y,gamma_temp] = EEGGenUnifNois(Gain_constrained,Vertices1,DeltaL,Ns,varsigrange,SrcExistRegion,SNR);%Generate EEG data contaminated by homoscedastic noise
% [Y,gamma_temp] = EEGGenNonunifNois(Gain_constrained,Vertices1,DeltaL,Ns,varsigrange,SrcExistRegion,SNR);%Generate EEG data contaminated by heteroscedastic noise
% [Y,gamma_temp] = EEGGenPartRegion(Gain_constrained,Vertices1,DeltaL,Ns,varsigrange,SrcExistRegion,SNR);%Generate EEG data (where half of voxels in ROIs are activated) contaminated by homoscedastic noise
% [Y,gamma_temp,VoxelGroundTruth] = EEGGenPartRegionNN(Gain_constrained,C,Vertices1,DeltaL,Ns,varsigrange,SrcExistRegion,SNR);%Generate EEG data (where half of voxels in ROIs are activated) contaminated by heteroscedastic noise
gamma(SrcExistRegion(:)) = gamma_temp;
% Y can be replaced by real EEG data to perform EEG cortical source localization algorithms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prposed: NNBSBL
[x1,VarVoxelNNBSBL] = NNBSBL(Y,Am_estimate,C,Vertices1);%VarVoxelNNBSBL consists of estimated voxel-level source variances.
AP1 = APrime(abs(x1),SrcExistRegion);
maxx1 = max(x1);
x1 = x1/maxx1;
dbx1 = 10*log(x1)/log(10);
maxdb1 = ceil(max(abs(dbx1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BSBL-2S
[x2,Xbsbl2s,VarVoxelBSBL2S] = BSBL2S(Y,U,s2,H,C,Gain_constrained,Vertices1);
AP2 = APrime(x2,SrcExistRegion);
maxx2 = max(x2);
x2 = x2/maxx2;
dbx2 = 10*log(x2)/log(10);
maxdb2 = ceil(max(abs(dbx2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Smooth Champagne with Robust Noise Estimation
[x3,XChampagneEM,VarVoxelSmoothChampagne] = ChampagneEM(Y,H,C,Hg,Gain_constrained,Vertices1);
AP3 = APrime(x3,SrcExistRegion);
maxx3 = max(x3);
x3 = x3/maxx3;
dbx3 = 10*log(x3)/log(10);
maxdb3 = ceil(max(abs(dbx3)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sLORETA
[x4,XsL,VarVoxelsLORETA] = sLORETA(Y,Gain_constrained,Vertices1,C);
AP4 = APrime(x4,SrcExistRegion);
maxx4 = max(x4);
x4 = x4/maxx4;
dbx4 = 10*log(x4)/log(10);
maxdb4 = ceil(max(abs(dbx4)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LCMV beamforming
[x5,XsLCMV,VarVoxelLCMV] = LCMVbeamforming(Y,Gain_constrained,Vertices1,C);
AP5 = APrime(x5,SrcExistRegion);
maxx5 = max(x5);
x5 = x5/maxx5;
dbx5 = 10*log(x5)/log(10);
maxdb5 = ceil(max(abs(dbx5)));

maxdb = max([maxdb1,maxdb2,maxdb3,maxdb4,maxdb5]);

figure(1);
subplot(5,1,1);
plot(1:K,dbx1,'rx-')
for k = 1:source_num
    hold on
    xv = SrcExistRegion(k)*ones(1,maxdb+1);
    yv = -maxdb:1:0;
    plot(xv,yv,'k:');
end
hold off
axis([0 150 -maxdb 0])
xlabel('ROI Index')
ylabel('Normalized var param (dB)')
legend('NNBSBL')
grid on
subplot(5,1,2);
plot(1:K,dbx2,'bx-')
for k = 1:source_num
    hold on
    xv = SrcExistRegion(k)*ones(1,maxdb+1);
    yv = -maxdb:1:0;
    plot(xv,yv,'k:');
end
hold off
axis([0 150 -maxdb 0])
xlabel('ROI Index')
ylabel('Normalized var param (dB)')
legend('BSBL-2S')
grid on
subplot(5,1,3);
plot(1:K,dbx3,'mx-')
for k = 1:source_num
    hold on
    xv = SrcExistRegion(k)*ones(1,maxdb+1);
    yv = -maxdb:1:0;
    plot(xv,yv,'k:');
end
hold off
axis([0 150 -maxdb 0])
xlabel('ROI Index')
ylabel('Normalized var param (dB)')
legend('Smooth Champagne')
grid on
subplot(5,1,4);
plot(1:K,dbx4,'gx-')
for k = 1:source_num
    hold on
    xv = SrcExistRegion(k)*ones(1,maxdb+1);
    yv = -maxdb:1:0;
    plot(xv,yv,'k:');
end
hold off
axis([0 150 -maxdb 0])
xlabel('ROI Index')
ylabel('Normalized var param (dB)')
legend('sLORETA')
grid on
subplot(5,1,5);
plot(1:K,dbx5,'kx-')
for k = 1:source_num
    hold on
    xv = SrcExistRegion(k)*ones(1,maxdb+1);
    yv = -maxdb:1:0;
    plot(xv,yv,'k:');
end
hold off
axis([0 150 -maxdb 0])
xlabel('ROI Index')
ylabel('Normalized var param (dB)')
legend('LCMV')
grid on

rmseloc1 = locerroreval(dbx1,SrcExistRegion,rangecore);
rmseloc2 = locerroreval(dbx2,SrcExistRegion,rangecore);
rmseloc3 = locerroreval(dbx3,SrcExistRegion,rangecore);
rmseloc4 = locerroreval(dbx4,SrcExistRegion,rangecore);
rmseloc5 = locerroreval(dbx5,SrcExistRegion,rangecore);

% save('SourceLoc.mat','VarVoxelNNBSBL','VarVoxelBSBL2S','VarVoxelSmoothChampagne','VarVoxelsLORETA','VarVoxelLCMV');