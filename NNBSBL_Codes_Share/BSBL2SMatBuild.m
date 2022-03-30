function [U,s2,H,Hg] = BSBL2SMatBuild(Gain_constrained,Vertices1,DeltaL)
Ng = length(DeltaL);
M = size(Gain_constrained,1);
H = cell(Ng,1);
A = [];
Hg = [];
for i = 1:Ng
    Vertices_temp = Vertices1{i};
    Htemp = Gain_constrained(:,Vertices_temp)*pinv(DeltaL{i});
    H{i} = Htemp;
    A = [A Htemp];
    Hg = [Hg Gain_constrained(:,Vertices_temp)];
end
[U, S, V] = svd(A);
s2 = (diag(S(1:M,1:M))).^2;