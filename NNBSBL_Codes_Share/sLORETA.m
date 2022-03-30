function [gamma,XsL,Wr] = sLORETA(Y,K,Vertices1,C)
alpha = 0.1;
[M,N] = size(K);
Ns = length(Y(1,:));
Ng = length(Vertices1);
gamma = zeros(Ng,1);
T = K.'*pinv(K*K.'+alpha*eye(M));
J = T*Y;
Sj = T*K;
XsL = zeros(N,Ns);
Wr = zeros(1,N);
for n = 1:N
    XsL(n,:) = J(n,:)/sqrt(Sj(n,n));
    Wr(n) = abs(XsL(n,:)*XsL(n,:).')/Ns;
end
Wr = Wr/max(Wr);

for i = 1:Ng
    Xtemp = XsL(Vertices1{i},:);
    Rx = Xtemp*Xtemp.'/Ns;
    gamma(i) = norm(Rx,'fro')/norm(C{i},'fro');
end
