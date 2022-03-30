clear all
close all
clc

ParamMat_name_data = dir([pwd,'\ParamMat\*.mat']);
ParamMat_dir_data = [pwd,'\ParamMat\'];
for kk = 1:size(ParamMat_name_data,1)
    ParamMatName = ParamMat_name_data(kk,1).name;
    filename = [ParamMat_dir_data ParamMatName];
    load(filename);
end
Gain_constrained([19,40],:)=[];
%---------------------------------------------------------------------------------%
[DeltaL,C] = LaplaceOperMatGen(Faces,GridLoc,Vertices1);
rangecore = brainrangecluster(GridLoc,Vertices1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Am_estimate=[];
regionweight = zeros(length(DeltaL),1);
for i=1:length(DeltaL)
    Vertices_temp=Vertices1{i};
    Atemp=Gain_constrained(:,Vertices_temp)*C{i}*Gain_constrained(:,Vertices_temp)';
    regionweight(i) = norm(Atemp,'fro');
    Atemp=Atemp/regionweight(i);
    A(i,:,:)=Atemp;
end
M = size(Gain_constrained,1);
for m=1:M
    A_m_temp=[];
    for i=1:length(DeltaL)
        A_temp=A(i,:,m);
        A_m_temp(:,i)= A_temp;
    end
    Am_estimate(:,:,m)=A_m_temp;
end
clear A A_m_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U,s2,H,Hg] = BSBL2SMatBuild(Gain_constrained,Vertices1,DeltaL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SampleRate = 1000;
SNR = 5;%SNR
varsigrange = [1,1];
source_num = 5;
K = length(Region);
MonteCarloNum = 1000;
numsamp = 60:10:120;%Considerd number of samples
lennumsamp = length(numsamp);
APsampNNBSBL = zeros(lennumsamp,1);
APsampBSBL2S = zeros(lennumsamp,1);
APsampChampagneEM = zeros(lennumsamp,1);
APsampsLORETA = zeros(lennumsamp,1);
APsampsLCMV = zeros(lennumsamp,1);
sampnum = 1;
for Ns = 60:10:120
    for mc = 1:MonteCarloNum
        Ns
        mc
        tic
        ind_temp = randperm(K,source_num);
        SrcExistRegion = sort(ind_temp,'ascend');
        gamma = zeros(K,1);
        [Y,gamma_temp] = EEGGenUnifNois(Gain_constrained,Vertices1,DeltaL,Ns,varsigrange,SrcExistRegion,SNR);%Generate EEG data contaminated by homoscedastic noise
        % [Y,gamma_temp] = EEGGenNonunifNois(Gain_constrained,Vertices1,DeltaL,Ns,varsigrange,SrcExistRegion,SNR);%Generate EEG data contaminated by heteroscedastic noise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %NNBSBL
        x1 = NNBSBL(Y,Am_estimate,C,Vertices1);
        AP1 = APrime(abs(x1),SrcExistRegion);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %BSBL-2S
        [x2,Xbsbl2s] = BSBL2S(Y,U,s2,H,C,Gain_constrained,Vertices1);
        AP2 = APrime(abs(x2),SrcExistRegion);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Smooth Champagne with Robust Noise Estimation
        [x3,XChampagneEM] = ChampagneEM(Y,H,C,Hg,Gain_constrained,Vertices1);
        AP3 = APrime(abs(x3),SrcExistRegion);
        %sLORETA
        [x4,XsL] = sLORETA(Y,Gain_constrained,Vertices1,C);
        AP4 = APrime(abs(x4),SrcExistRegion);
        %LCMV
        [x5,XsLCMV] = LCMVbeamforming(Y,Gain_constrained,Vertices1,C);
        AP5 = APrime(abs(x5),SrcExistRegion);
        
        APsampNNBSBL(sampnum) = APsampNNBSBL(sampnum) + AP1;
        APsampBSBL2S(sampnum) = APsampBSBL2S(sampnum) + AP2;
        APsampChampagneEM(sampnum) = APsampChampagneEM(sampnum) + AP3;
        APsampsLORETA(sampnum) = APsampsLORETA(sampnum) + AP4;
        APsampsLCMV(sampnum) = APsampsLCMV(sampnum) + AP5;
        toc
    end
    APsampNNBSBL(sampnum) = APsampNNBSBL(sampnum)/MonteCarloNum;
    APsampBSBL2S(sampnum) = APsampBSBL2S(sampnum)/MonteCarloNum;
    APsampChampagneEM(sampnum) = APsampChampagneEM(sampnum)/MonteCarloNum;
    APsampsLORETA(sampnum) = APsampsLORETA(sampnum)/MonteCarloNum;
    APsampsLCMV(sampnum) = APsampsLCMV(sampnum)/MonteCarloNum;
    sampnum = sampnum + 1;
end
figure(1);
plot(numsamp,APsampNNBSBL,'ro-');
hold on
plot(numsamp,APsampBSBL2S,'bx-');
hold on
plot(numsamp,APsampChampagneEM,'m>-');
hold on
plot(numsamp,APsampsLORETA,'ks-');
hold on
plot(numsamp,APsampsLCMV,'yd-');
hold off
axis([numsamp(1) numsamp(end) 0 1])
xlabel('Number of Samples');
ylabel('A Prime');
legend('NNBSBL','BSBL-2S','Smooth Champagne','sLORETA','LCMV');
set(gca,'XTick',60:10:120)  
set(gca,'XTickLabel',{'60','70','80','90','100','110','120'})