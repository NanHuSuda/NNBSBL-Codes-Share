function [gamma,XsL,Wr] = LCMVbeamforming(Y,K,Vertices1,C)
alpha = 1e-5;
[M,N] = size(K);
Ns = length(Y(1,:));
Ng = length(Vertices1);
gamma = zeros(Ng,1);
pCov = pinv(Y*Y.'/Ns+alpha*eye(M));
Atemp = zeros(M,1);
XsL = zeros(N,Ns);
Wr = zeros(1,N);
for n = 1:N
    Atemp(:) = K(:,n);
    XsL(n,:) = Atemp.'*pCov*Y/(Atemp.'*pCov*Atemp);
    Wr(n) = abs(XsL(n,:)*XsL(n,:).')/Ns;
end
Wr = Wr/max(Wr);

for i = 1:Ng
    Xtemp = XsL(Vertices1{i},:);
    Rx = Xtemp*Xtemp.'/Ns;
    gamma(i) = norm(Rx,'fro')/norm(C{i},'fro');
end