function [DeltaL,C,L] = LaplaceOperMatGen(Faces,GridLoc,Vertices1)
%Faces��������ʽ��ÿ�б�ʾһ�������������������ģ��ڵ��̳������ж�Ӧ�еģ�����
%GridLoc��������ʽ����ʾΪ���̳�����ÿ�ж�Ӧ���������ֵ��ά��Ϊ ��������3
%Vertices1��cell��ʽ����Ԫ�ر�ʾ�����������������ж���ģ��ڵ��̳������ж�Ӧ�еģ��������ϣ���������ʽ��
Nvert = size(GridLoc,1);%ȫ������
Ntria = size(Faces,1);%�����������εĸ���
H = zeros(Nvert,Nvert);%���������ϵ����ڵ��ľ��룻0��ʾ�����ڣ�����˳��Ϊ���̳�������е�����˳��
%�������ڵ�������H
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
%�����ڵ�������H����ȫ����������˹���Ӿ���L����ά��Ϊ��ȫ��������ȫ������
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
%��ȫ����������˹���Ӿ���L����ȡһ��cell�ṹ����C������ÿ��Ԫ����Ԫ�洢����ÿ��������Laplace���Ӿ���
%Ԫ����Ԫ������˳����Region��������˳��Ԫ����Ԫ��ÿ����������������˳����Vertices1����˳��
Nregi = size(Vertices1,1);%��������
DeltaL = cell(Nregi,1);
C = cell(Nregi,1);
for n = 1:Nregi
    DeltaL{n} = L(Vertices1{n},Vertices1{n});
    C{n} = pinv(DeltaL{n}.'*DeltaL{n});
end