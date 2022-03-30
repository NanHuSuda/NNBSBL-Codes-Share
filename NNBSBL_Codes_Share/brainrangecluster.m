function rangecore = brainrangecluster(GridLoc,Vertices1)
Ng = length(Vertices1);
rangecore = zeros(3,Ng);
for i = 1:Ng
    Ni = length(Vertices1{i});
    rangecore(:,i) = (2/Ni)*mean(GridLoc(Vertices1{i},:));
end
    