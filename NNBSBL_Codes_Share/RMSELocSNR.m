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
dattimlen = 0.080;%data length used
Ns = floor(SampleRate*dattimlen);
varsigrange = [1,1];
source_num = 3;
K = length(Region);

MonteCarloNum = 1;
snr = -10:3:20;
lensnr = length(snr);
rmsesnrNNBSBL = zeros(lensnr,1);
rmsesnrBSBL2S = zeros(lensnr,1);
rmsesnrChampagneEM = zeros(lensnr,1);
rmsesnrsLORETA = zeros(lensnr,1);
rmsesnrLCMV = zeros(lensnr,1);
snrnum = 1;
for SNR = -10:3:20
    for mc = 1:MonteCarloNum
        SNR
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
        %¶ÔBSBL-2S
        [x2,Xbsbl2s] = BSBL2S(Y,U,s2,H,C,Gain_constrained,Vertices1);
        maxx2 = max(x2);
        x2 = x2/maxx2;
        dbx2 = 10*log(x2)/log(10);
        maxdb2 = ceil(max(abs(dbx2)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Smooth Champagne
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
        rmsesnrNNBSBL(snrnum) = rmsesnrNNBSBL(snrnum) + rmseloc1^2;
        rmsesnrBSBL2S(snrnum) = rmsesnrBSBL2S(snrnum) + rmseloc2^2;
        rmsesnrChampagneEM(snrnum) = rmsesnrChampagneEM(snrnum) + rmseloc3^2;
        rmsesnrsLORETA(snrnum) = rmsesnrsLORETA(snrnum) + rmseloc4^2;
        rmsesnrLCMV(snrnum) = rmsesnrLCMV(snrnum) + rmseloc5^2;
        toc
    end
    rmsesnrNNBSBL(snrnum) = sqrt(rmsesnrNNBSBL(snrnum)/MonteCarloNum);
    rmsesnrBSBL2S(snrnum) = sqrt(rmsesnrBSBL2S(snrnum)/MonteCarloNum);
    rmsesnrChampagneEM(snrnum) = sqrt(rmsesnrChampagneEM(snrnum)/MonteCarloNum);
    rmsesnrsLORETA(snrnum) = sqrt(rmsesnrsLORETA(snrnum)/MonteCarloNum);
    rmsesnrLCMV(snrnum) = sqrt(rmsesnrLCMV(snrnum)/MonteCarloNum);
    snrnum = snrnum + 1;
end
rmsesnrNNBSBL = 100*rmsesnrNNBSBL;
rmsesnrBSBL2S = 100*rmsesnrBSBL2S;
rmsesnrChampagneEM = 100*rmsesnrChampagneEM;
rmsesnrsLORETA = 100*rmsesnrsLORETA;
rmsesnrLCMV = 100*rmsesnrLCMV;

figure(1);
plot(snr,rmsesnrNNBSBL,'ro-');
hold on
plot(snr,rmsesnrBSBL2S,'bx-');
hold on
plot(snr,rmsesnrChampagneEM,'m>-');
hold on
plot(snr,rmsesnrsLORETA,'ks-');
hold on
plot(snr,rmsesnrLCMV,'yd-');
hold off
axis([snr(1) snr(end) 0 inf])
xlabel('SNR (dB)');
ylabel('RMSE of region locations (cm)');
legend('NNBSBL','BSBL-2S','Smooth Champagne','sLORETA','LCMV');
set(gca,'XTick',-10:3:20)  
set(gca,'XTickLabel',{'-10','-7','-4','-1','2','5','8','11','14','17','20'})