function [gamma,X,VarVoxel] = ChampagneEM(Y,H,C,Hg,Gain_constrained,Vertices1)
epsilon = 0.1;
[Ny,N] = size(Y);
Ng = length(H);
Cy = Y*Y.'/N;
%initialization
GammaNois_init = 0.01*norm(Cy,'fro')/Ny*ones(Ny,1);%noise covariance matrix
F = [];
for i = 1:Ng
    F = [F H{i}];
end
X0 = F.'*pinv(F*F.'+1e-5*eye(Ny))*Y;
gamma_init = zeros(Ng,1);
NGrid = 0;
for i = 1:Ng
    Ni = size(H{i},2);
    NGrid = NGrid + Ni;
    gamma_init(i) = (norm(X0(1:Ni,:),'fro'))^2/(Ni*N);
    X0(1:Ni,:) = [];
end
%iteration
GammaNois = GammaNois_init;
gamma = gamma_init;
Sigmay = diag(GammaNois);
for i = 1:Ng
    Sigmay = Sigmay + gamma(i)*H{i}*H{i}.';
end
Sigmayinv = pinv(Sigmay);
modelevidence = trace(Cy*Sigmayinv);
% modelevidence = log(abs(det(Sigmay))) + trace(Cy*Sigmayinv);%model evidence by initialization
modelevidence_old = -1000;
NumChampIter = 0;
while abs(modelevidence-modelevidence_old) >= epsilon
    NumChampIter = NumChampIter + 1;
    modelevidence_old = modelevidence;
    gamma_old = gamma;
    %calculate mean and covariance of signals at voxels: X and Cx
    SigmaxHT = [];
    SigmaxG = zeros(NGrid,NGrid);
    numbegin = 1;
    for i = 1:Ng
        Ni = size(H{i},2);
        SigmaxG(numbegin:numbegin+Ni-1,numbegin:numbegin+Ni-1) = gamma(i)*C{i};
        numbegin = numbegin+Ni;
        SigmaxHT = [SigmaxHT;gamma(i)*C{i}*Gain_constrained(:,Vertices1{i}).'];
    end
    X = SigmaxHT*Sigmayinv*Y;%mean matrix
    Sigmax = SigmaxG - SigmaxHT*Sigmayinv*SigmaxHT.';
    %update variances in the noise covariance matrix
    for k = 1:Ny
        GammaNois(k) = 0;
        for n = 1:N
            GammaNois(k) = GammaNois(k) + (1/N)*(Y(k,n)-Hg(k,:)*X(:,n))^2;
        end
        GammaNois(k) = GammaNois(k) + Hg(k,:)*Sigmax*Hg(k,:).';
    end
    %update variances in each brain region
    for i = 1:Ng
        gamma(i) = (gamma_old(i)/sqrt(N))*norm(H{i}.'*Sigmayinv*Y,'fro')/sqrt(trace(H{i}.'*Sigmayinv*H{i}));
    end
    Sigmay = diag(GammaNois);
    for i = 1:Ng
        Sigmay = Sigmay + gamma(i)*H{i}*H{i}.';
    end
    Sigmayinv = pinv(Sigmay);
%     modelevidence = log(abs(det(Sigmay))) + trace(Cy*Sigmayinv)
    modelevidence = trace(Cy*Sigmayinv);
end
%Source estimation
SigmaxHT = [];
SigmaxG = zeros(NGrid,NGrid);
numbegin = 1;
for i = 1:Ng
    Ni = size(H{i},2);
    SigmaxG(numbegin:numbegin+Ni-1,numbegin:numbegin+Ni-1) = gamma(i)*C{i};
    numbegin = numbegin+Ni;
    SigmaxHT = [SigmaxHT;gamma(i)*C{i}*Gain_constrained(:,Vertices1{i}).'];
end
X = SigmaxHT*Sigmayinv*Y;%mean matrix
Sigmax = SigmaxG - SigmaxHT*Sigmayinv*SigmaxHT.';

ss = gamma;
weight=[];
ind=[];
for i=1:Ng 
    weight_temp=diag(C{i})*ss(i);
    weight=[weight;weight_temp];
    ind_temp=Vertices1{i};
    ind=[ind,ind_temp];
end

weight_new(ind)=weight;
VarVoxel=weight_new/max(weight_new);
