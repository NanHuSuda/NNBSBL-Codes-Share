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
SNR = 5;
dattimlen = 0.080;
Ns = floor(SampleRate*dattimlen);
varsigrange = [1,1];
K = length(Region);
MonteCarloNum = 1;
SorNum = 2:8;
lensor = length(SorNum);
APsnrNNBSBL = zeros(lensor,1);
APsesnrBSBL2S = zeros(lensor,1);
APsnrChampagneEM = zeros(lensor,1);
APsnrsLORETA = zeros(lensor,1);
APsorLCMV = zeros(lensor,1);
sornum = 1;
for source_num = 2:8
    for mc = 1:MonteCarloNum
        source_num
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %sLORETA
        [x4,XsL] = sLORETA(Y,Gain_constrained,Vertices1,C);
        AP4 = APrime(abs(x4),SrcExistRegion);
        %LCMV
        [x5,XsLCMV] = LCMVbeamforming(Y,Gain_constrained,Vertices1,C);
        AP5 = APrime(abs(x5),SrcExistRegion);
        APsnrNNBSBL(sornum) = APsnrNNBSBL(sornum) + AP1;
        APsesnrBSBL2S(sornum) = APsesnrBSBL2S(sornum) + AP2;
        APsnrChampagneEM(sornum) = APsnrChampagneEM(sornum) + AP3;
        APsnrsLORETA(sornum) = APsnrsLORETA(sornum) + AP4;
        APsorLCMV(sornum) = APsorLCMV(sornum) + AP5;
        toc
    end
    APsnrNNBSBL(sornum) = APsnrNNBSBL(sornum)/MonteCarloNum;
    APsesnrBSBL2S(sornum) = APsesnrBSBL2S(sornum)/MonteCarloNum;
    APsnrChampagneEM(sornum) = APsnrChampagneEM(sornum)/MonteCarloNum;
    APsnrsLORETA(sornum) = APsnrsLORETA(sornum)/MonteCarloNum;
    APsorLCMV(sornum) =APsorLCMV(sornum)/MonteCarloNum;
    sornum = sornum + 1;
end

figure(1);
plot(SorNum,APsnrNNBSBL,'ro-');
hold on
plot(SorNum,APsesnrBSBL2S,'bx-');
hold on
plot(SorNum,APsnrChampagneEM,'m>-');
hold on
plot(SorNum,APsnrsLORETA,'ks-');
hold on
plot(SorNum,APsorLCMV,'yd-');
hold off
axis([SorNum(1) SorNum(end) 0 1])
xlabel('Number of Source Existing Regions');
ylabel('A Prime');
legend('NNBSBL','BSBL-2S','Smooth Champagne','sLORETA','LCMV');
set(gca,'XTick',2:8)  
set(gca,'XTickLabel',{'2','3','4','5','6','7','8'})