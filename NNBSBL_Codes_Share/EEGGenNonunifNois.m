function [Y,gamma] = EEGGenNonunifNois(Gain_constrained,Vertices1,DeltaL,Ns,varsigrange,SrcExistRegion,SNR)
K = length(SrcExistRegion);
gamma = varsigrange(1)*ones(K,1)+(varsigrange(2)-varsigrange(1))*rand(K,1);
A = [];
S = [];
for k = 1:K
    Vertices_temp = Vertices1{SrcExistRegion(k)};
    A = [A Gain_constrained(:,Vertices_temp(:))];
    DeltaL_temp = DeltaL{SrcExistRegion(k)};
    S_temp = sqrt(gamma(k))*pinv(DeltaL_temp)*randn(length(Vertices_temp),Ns);
    S = [S;S_temp];
end
Ynf = A*S;
randnoisvar = 0.5*ones(size(Gain_constrained,1),1)+1*rand(size(Gain_constrained,1),1);
Nmat = randn(size(Gain_constrained,1),Ns).*(randnoisvar*ones(1,Ns));
lamda = 10^(-SNR/10)*(norm(Ynf,'fro')/norm(Nmat,'fro'))^2;
Y = Ynf + sqrt(lamda)*Nmat;