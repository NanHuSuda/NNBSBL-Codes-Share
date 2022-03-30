function AP = APrime(x,SrcExistRegion)
Thamp = 0.01;%a ratio multiplying the maximum value to give a threshold
Thratio = 0.1;%the chosen ratio of data in descending order
Ng = length(x);%total number of ROIs
K = length(SrcExistRegion);%number of activated ROIs
xmval = max(x);
[xd,regd] = dsort(x);
detregcollection = [];
for k = 1:ceil(Thratio*Ng)
    if xd(k) >= Thamp*xmval
        detregcollection = [detregcollection;regd(k)];
    end
end
Nhit = 0;
Nd = length(detregcollection);
for k = 1:K
    for n = 1:Nd
        if detregcollection(n) == SrcExistRegion(k)
            Nhit = Nhit + 1;
        end
    end
end
Nfal = Nd - Nhit;
Hr = Nhit/K;
Fr = Nfal/(Ng-K);
AP = (Hr-Fr)/2+0.5;
