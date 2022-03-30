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
SNR = 5;%SNR
varsigrange = [1,1];
source_num = 3;
K = length(Region);
MonteCarloNum = 1000;
numsamp = 60:10:120;
lennumsamp = length(numsamp);
rmsesampNNBSBL = zeros(lennumsamp,1);
rmsesampBSBL2S = zeros(lennumsamp,1);
rmsesampChampagneEM = zeros(lennumsamp,1);
rmsesampsLORETA = zeros(lennumsamp,1);
rmsesampLCMV = zeros(lennumsamp,1);
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
        maxx1 = max(x1);
        x1 = x1/maxx1;
        dbx1 = 10*log(x1)/log(10);
        maxdb1 = ceil(max(abs(dbx1)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %对比算法：BSBL-2S
        [x2,Xbsbl2s] = BSBL2S(Y,U,s2,H,C,Gain_constrained,Vertices1);
        maxx2 = max(x2);
        x2 = x2/maxx2;
        dbx2 = 10*log(x2)/log(10);
        maxdb2 = ceil(max(abs(dbx2)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %对比算法：BSBL-2S
        [x3,XChampagneEM] = ChampagneEM(Y,H,C,Hg,Gain_constrained,Vertices1);
        maxx3 = max(x3);
        x3 = x3/maxx3;
        dbx3 = 10*log(x3)/log(10);
        maxdb3 = ceil(max(abs(dbx3)));
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %sLORETA
        [x4,XsL] = sLORETA(Y,Gain_constrained,Vertices1,C);
        maxx4 = max(x4);
        x4 = x4/maxx4;
        dbx4 = 10*log(x4)/log(10);
        maxdb4 = ceil(max(abs(dbx4)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %LCMV
        [x5,XsLCMV] = LCMVbeamforming(Y,Gain_constrained,Vertices1,C);
        maxx5 = max(x5);
        x5 = x5/maxx5;
        dbx5 = 10*log(x5)/log(10);
        maxdb5 = ceil(max(abs(dbx5)));
        
        maxdb = max([maxdb1,maxdb2,maxdb3,maxdb4,maxdb5]);

        dbx1 = dbx1*(maxdb/maxdb1);
        dbx2 = dbx2*(maxdb/maxdb2);
        dbx3 = dbx3*(maxdb/maxdb3);
        dbx4 = dbx4*(maxdb/maxdb4);
        dbx5 = dbx5*(maxdb/maxdb5);
        rmseloc1 = locerroreval(dbx1,SrcExistRegion,rangecore);
        rmseloc2 = locerroreval(dbx2,SrcExistRegion,rangecore);
        rmseloc3 = locerroreval(dbx3,SrcExistRegion,rangecore);
        rmseloc4 = locerroreval(dbx4,SrcExistRegion,rangecore);
        rmseloc5 = locerroreval(dbx5,SrcExistRegion,rangecore);
        rmsesampNNBSBL(sampnum) = rmsesampNNBSBL(sampnum) + rmseloc1^2;
        rmsesampBSBL2S(sampnum) = rmsesampBSBL2S(sampnum) + rmseloc2^2;
        rmsesampChampagneEM(sampnum) = rmsesampChampagneEM(sampnum) + rmseloc3^2;
        rmsesampsLORETA(sampnum) = rmsesampsLORETA(sampnum) + rmseloc4^2;
        rmsesampLCMV(sampnum) = rmsesampLCMV(sampnum) + rmseloc5^2;
        toc
    end
    rmsesampNNBSBL(sampnum) = sqrt(rmsesampNNBSBL(sampnum)/MonteCarloNum);
    rmsesampBSBL2S(sampnum) = sqrt(rmsesampBSBL2S(sampnum)/MonteCarloNum);
    rmsesampChampagneEM(sampnum) = sqrt(rmsesampChampagneEM(sampnum)/MonteCarloNum);
    rmsesampsLORETA(sampnum) = sqrt(rmsesampsLORETA(sampnum)/MonteCarloNum);
    rmsesampLCMV(sampnum) = sqrt(rmsesampLCMV(sampnum)/MonteCarloNum);
    sampnum = sampnum + 1;
end
figure(1);
rmsesampNNBSBL = 100*rmsesampNNBSBL;
rmsesampBSBL2S = 100*rmsesampBSBL2S;
rmsesampChampagneEM = 100*rmsesampChampagneEM;
rmsesampsLORETA = 100*rmsesampsLORETA;
rmsesampLCMV = 100*rmsesampLCMV;
plot(numsamp,rmsesampNNBSBL,'ro-');
hold on
plot(numsamp,rmsesampBSBL2S,'bx-');
hold on
plot(numsamp,rmsesampChampagneEM,'m>-');
hold on
plot(numsamp,rmsesampsLORETA,'ks-');
hold on
plot(numsamp,rmsesampLCMV,'yd-');
hold off
axis([numsamp(1) numsamp(end) 0 inf])
xlabel('Number of Samples');
ylabel('RMSE of region locations');
legend('NNBSBL','BSBL-2S','Smooth Champagne','sLORETA','LCMV');
set(gca,'XTick',60:10:120)  
set(gca,'XTickLabel',{'60','70','80','90','100','110','120'})