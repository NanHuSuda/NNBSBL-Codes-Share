function [DeltaL,C,L] = LaplaceOperMatGen(Faces,GridLoc,Vertices1)
%Faces：矩阵形式，每行表示一个三角网格的三个顶点的（在导程场矩阵中对应列的）索引
%GridLoc：矩阵形式，表示为导程场矩阵每列对应顶点的坐标值，维度为 顶点数×3
%Vertices1：cell形式，各元素表示该脑区所包含的所有顶点的（在导程场矩阵中对应列的）索引集合（行向量形式）
Nvert = size(GridLoc,1);%全顶点数
Ntria = size(Faces,1);%网格中三角形的个数
H = zeros(Nvert,Nvert);%三角网格上的相邻点间的距离；0表示不相邻；索引顺序为导程场矩阵各列的索引顺序
%计算相邻点距离矩阵H
for k = 1:Ntria
    if H(Faces(k,1),Faces(k,2)) == 0
        H(Faces(k,1),Faces(k,2)) = norm(GridLoc(Faces(k,1))-GridLoc(Faces(k,2)),2);
        H(Faces(k,2),Faces(k,1)) = H(Faces(k,1),Faces(k,2));
    end
    if H(Faces(k,1),Faces(k,3)) == 0
        H(Faces(k,1),Faces(k,3)) = norm(GridLoc(Faces(k,1))-GridLoc(Faces(k,3)),2);
        H(Faces(k,3),Faces(k,1)) = H(Faces(k,1),Faces(k,3));
    end
    if H(Faces(k,2),Faces(k,3)) == 0
        H(Faces(k,2),Faces(k,3)) = norm(GridLoc(Faces(k,2))-GridLoc(Faces(k,3)),2);
        H(Faces(k,3),Faces(k,2)) = H(Faces(k,2),Faces(k,3));
    end
end
%由相邻点距离矩阵H计算全顶点拉普拉斯算子矩阵L，其维度为：全顶点数×全顶点数
hm = 1./sum(H)';
Hr = H;
for k1 = 1:Nvert
    for k2 = 1:Nvert
        if Hr(k1,k2)
            Hr(k1,k2) = 1/Hr(k1,k2);
        end
    end
end
hr = sum(Hr)';
L = 4*Hr.*repmat(hm,1,Nvert);
for k = 1:Nvert
    L(k,k) = -4*hm(k)*hr(k);
end
%从全顶点拉普拉斯算子矩阵L中提取一个cell结构数据C，其中每个元胞单元存储的是每个脑区的Laplace算子矩阵
%元胞单元的排列顺序按照Region中脑区的顺序，元胞单元内每个矩阵行与列排列顺序按照Vertices1排列顺序
Nregi = size(Vertices1,1);%脑区个数
DeltaL = cell(Nregi,1);
C = cell(Nregi,1);
for n = 1:Nregi
    DeltaL{n} = L(Vertices1{n},Vertices1{n});
    C{n} = pinv(DeltaL{n}.'*DeltaL{n});
end