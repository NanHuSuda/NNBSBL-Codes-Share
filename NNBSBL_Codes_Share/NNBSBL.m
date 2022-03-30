function [ss,VarVoxel] = NNBSBL(X,A,C,Vertices1)
X=X/max(max(abs(X)));
j = sqrt(-1);
[M,T] = size(X);
N = size(A,2);
R = X*X'/T;%Sample covariance matrix
Yb = zeros(M-1,M);
Fyb = zeros(M-1,N,M);
for m = 1:M
    ytemp = R(:,m);
    ytemp(m) = [];
    Fytemp = zeros(M,N);
    Fytemp(:,:) = A(:,:,m);
    Fytemp(m,:) = [];
    Rwm = R(m,m)*R/T;
    Rwm(m,:) = [];
    Rwm(:,m) = [];
    Wm = inv(chol(Rwm+1e-6*eye(M-1),'lower'));
    Yb(:,m) = Wm*ytemp(:);
    Fyb(:,:,m) = Wm*Fytemp;
end
%%
iternum = 1000;%the number of iterations
L = M;
M2 = M-1;
% Initialization
z = zeros(N,1);%transformed data
H = zeros(N,N);
gamma = zeros(N,1);
for k = 1:L
    Fy = zeros(M2,N);
    Fy(:,:) = Fyb(:,:,k);
    y = zeros(M2,1);
    y(:) = Yb(:,k);
    z = z + Fy.'*y;
    H = H + Fy.'*Fy;
    gamma = gamma + abs(y'*Fy*Fy.'*y)/norm(Fy.'*Fy*Fy.'*Fy,'fro')*ones(N,1);
end
gamma0 = gamma/L;

multparam = [0.0005, 0.002, 0.01, 0.05, 0.1, 1, 10, 50, 200, 2000];
Np = length(multparam);
Qval = zeros(Np,1);
gq = zeros(N,Np);
for t = 1:Np
    gamma = multparam(t)*gamma0;%Multiple Initialization. However, only one initialization "gamma0" is also enough to get a fast and accurate iterative algorithm.
    Gamma = diag(gamma);
    Gamma_old = 0.001*ones(size(A,2), size(A,2));
    %Proceed
    k=1;
    while (norm(Gamma-Gamma_old,2)/norm(Gamma_old,2) > 1e-4) && (k < iternum)
        Gamma_old = Gamma;
        Zigma = pinv(H+pinv(Gamma_old));
        mu = Zigma*z;
        for n = 1:N
            Gamma(n,n) = 0;
            if (-1*mu(n)/sqrt(2*abs(Zigma(n,n)))) < 10
                Gamma(n,n) = abs(Zigma(n,n)) + mu(n)^2 + mu(n)*sqrt(2*abs(Zigma(n,n))/pi)*exp(-1*mu(n)^2/(2*abs(Zigma(n,n))))/erfc(real(-1*mu(n)/sqrt(2*abs(Zigma(n,n)))));
            else
                Gamma(n,n) = abs(Zigma(n,n));
            end
        end
        k=k+1;
    end
    Zigma = pinv(H+pinv(Gamma));
    mu = Zigma*z;
    E1 = zeros(N,1);
    E2 = zeros(N,N);
    for n = 1:N
%         h = real(-1*mu(n)/sqrt(2*abs(Zigma(n,n))));
%         if h < 10
%             E1(n) = mu(n) + sqrt(2*Zigma(n,n)/pi)*exp(-1*mu(n)^2/(2*abs(Zigma(n,n))))/erfc(h);
%             E2(n,n) = abs(Zigma(n,n)) + mu(n)^2 + mu(n)*sqrt(2*abs(Zigma(n,n))/pi)*exp(-1*mu(n)^2/(2*abs(Zigma(n,n))))/erfc(h);
%         else
%             E1(n) = 0;
%             E2(n,n) = abs(Zigma(n,n));
%         end
        E1(n) = mu(n);
        E2(n,n) = abs(Zigma(n,n))+mu(n)*mu(n);
    end
    for m = 2:N
        for n = 1:m-1
            E2(m,n) = Zigma(m,n)+mu(m)*mu(n);
            E2(n,m) = Zigma(m,n)+mu(m)*mu(n);
        end
    end
%     for m = 2:N
%         for n = 1:m-1
%             mum = mu(m);
%             mun = mu(n);
%             zm2 = Zigma(m,m);
%             zn2 = Zigma(n,n);
%             zmn = Zigma(m,n);
%             try 
%                 cmn = mvncdf([mum mun],[0,0],[zm2 zmn;zmn zn2]);
%                 kmm = zn2/(zm2*zn2-zmn^2);
%                 knn = zm2/(zm2*zn2-zmn^2);
%                 E2(m,n) = zmn + mum*mun + ((mum*zn2+mun*zmn)*sqrt(pi/(2*kmm))*exp(-mun^2/(2*zm2))*erfc(-sqrt(kmm/2)*(mum-mun*zmn/zn2))+(mum*zmn+mun*zm2)*sqrt(pi/(2*knn))*exp(-mum^2/(2*zn2))*erfc(-sqrt(knn/2)*(mun-mum*zmn/zm2)))/(2*pi*sqrt(zm2*zn2-zmn^2)*cmn);
%                 E2(n,m) = E2(m,n);
%             catch
%                 E2(m,n,l) = zmn;
%                 E2(n,m,l) = zmn;
%             end
%         end
%     end
    Qval(t) = E1.'*z - 0.5*trace((H+pinv(Gamma))*E2) - 0.5*sum(log(diag(Gamma)));
    gq(:,t) = abs(diag(Gamma));
end
Qmaxidx = 0;
thesholdres1 = 1;
thesholdres2 = 1.5;
for k = 2:Np-1
    if (Qval(k)-Qval(k-1)>thesholdres1)&&(abs(Qval(k)-Qval(k+1))<thesholdres2)
        Qmaxidx = k;
    elseif (abs(Qval(k)-Qval(k-1))<thesholdres2)&&(Qval(k)-Qval(k+1)>thesholdres1)
        Qmaxidx = k;
        break;
    end
end
if Qmaxidx == 0
    [Qmax, Qmaxidx] = max(Qval);
end
ss = gq(:,Qmaxidx);
Zigma = pinv(H+pinv(diag(ss)));
mu = Zigma*z;
ss = zeros(N,1);
for n = 1:N
    h = real(-1*mu(n)/sqrt(2*abs(Zigma(n,n))));
    if h < 10
        ss(n) = mu(n) + sqrt(2*Zigma(n,n)/pi)*exp(-1*mu(n)^2/(2*abs(Zigma(n,n))))/erfc(h);
    else
        ss(n) = 1e-9;
    end
end 
ss = abs(ss);%Estimated variance parameters for ROIs

weight=[];
ind=[];
for i=1:size(A,2) 
    weight_temp=diag(C{i})*ss(i);
    weight=[weight;weight_temp];
    ind_temp=Vertices1{i};
    ind=[ind,ind_temp];
end

weight_new(ind)=weight;
VarVoxel=weight_new/max(weight_new);%Estimated voxel-level source variances

