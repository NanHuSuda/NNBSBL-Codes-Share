function [gamma,X,VarVoxel] = BSBL2S(Y,U,s2,H,C,Gain_constrained,Vertices1)
epsilon = 0.1;
[Ny,N] = size(Y);
Ng = length(H);
Cy = Y*Y.'/N;
%initialization
Yb = (U.'*Y).';
yb2 = zeros(Ny,1);
for i = 1:Ny
    yb2(i) = (1/N)*(norm(Yb(:,i),2))^2;
end
F = [s2 ones(Ny,1)];
parminit = inv(F.'*F)*F.'*yb2;
lamda = parminit(2);
gammaF = parminit(1);
%stage-1: LORETA
phi = gammaF*s2+lamda*ones(Ny,1);
modelevidence = sum(log(phi)) + sum(yb2.*phi.^(-1));
modelevidence_old = -1000;
NumStage1 = 0;
while abs(modelevidence-modelevidence_old) >= epsilon
    NumStage1 = NumStage1 + 1;
    modelevidence_old = modelevidence;
    lamda = lamda*sum(yb2.*phi.^(-2))/sum(phi.^(-1));
    gammaF = gammaF*sum(yb2.*s2.*phi.^(-2))/sum(s2.*phi.^(-1));
    phi = gammaF*s2+lamda*ones(Ny,1);
    modelevidence = sum(log(phi)) + sum(yb2.*phi.^(-1));
end
%stage-2: pruning
modelevidence_old = -1000;
gamma = gammaF*ones(Ng,1);
Sigmay = lamda*eye(Ny);
for i = 1:Ng
    Sigmay = Sigmay + gamma(i)*H{i}*H{i}.';
end
Sigmayinv = pinv(Sigmay);
NumStage2 = 0;
while abs(modelevidence-modelevidence_old) >= epsilon
    NumStage2 = NumStage2 + 1;
    modelevidence_old = modelevidence;
    for i = 1:Ng
        gamma(i) = (gamma(i)/sqrt(N))*norm(H{i}.'*Sigmayinv*Y,'fro')/sqrt(trace(H{i}.'*Sigmayinv*H{i}));
    end
    Sigmay = lamda*eye(Ny);
    for i = 1:Ng
        Sigmay = Sigmay + gamma(i)*H{i}*H{i}.';
    end
    Sigmayinv = pinv(Sigmay);
%     modelevidence = log(abs(det(Sigmay))) + trace(Cy*Sigmayinv);
    modelevidence = trace(Cy*Sigmayinv);
end
%Source estimation
SigmaxHT = [];
for i = 1:Ng
    SigmaxHT = [SigmaxHT;gamma(i)*C{i}*Gain_constrained(:,Vertices1{i}).'];
end
X = SigmaxHT*Sigmayinv*Y;

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

    