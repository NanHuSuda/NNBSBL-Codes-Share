function rmseloc = locerroreval(x,SrcExistRegion,rangecore)
%x: normalized variance parameter vector (dB)
xthreshold = -18;
K = length(SrcExistRegion);
Ng = length(rangecore(1,:));
peakcollection = [];
peakamp = [];
rmseloc = 0;
for n = 1:Ng
    if x(n) >= xthreshold
        peakcollection = [peakcollection;n];
        peakamp = [peakamp;x(n)];
    end
end
peakcollection = sort(peakcollection);
if ~isempty(peakcollection)
    if length(peakcollection) >= K
        if length(peakcollection) > 3*K
            [peakampd, dnum] = dsort(peakamp);
            peakcollection = peakcollection(dnum(1:3*K));
        end
        for k = 1:K
            tempdist = [];
            for i = 1:length(peakcollection)
                tempdist(i,1) = i;
                tempdist(i,2) = norm(rangecore(peakcollection(i))-rangecore(SrcExistRegion(k)),2)^2;
            end
            [distmin,nummin] = min(tempdist(:,2));
            rmseloc = rmseloc + distmin;
            peakcollection(tempdist(nummin,1)) = [];
        end
        rmseloc = sqrt(rmseloc/K);
    else
        %Find locations of multiple maximas
        peakcollection = [];
        [xd,numd] = sort(abs(x));
        peakcollection = numd(1:2*K);
        peakcollection = sort(peakcollection);
        for k = 1:K
            tempdist = [];
            for i = 1:length(peakcollection)
                tempdist(i,1) = i;
                tempdist(i,2) = norm(rangecore(peakcollection(i))-rangecore(SrcExistRegion(k)),2)^2;
            end
            [distmin,nummin] = min(tempdist(:,2));
            rmseloc = rmseloc + distmin;
            peakcollection(tempdist(nummin,1)) = [];
        end
        rmseloc = sqrt(rmseloc/K);
    end
end
            