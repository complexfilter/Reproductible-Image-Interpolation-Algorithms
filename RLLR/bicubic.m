function Znew=bicubic(f)
%该函数是正对数字图像处理中的双三次线性插值函数
%f是所要处理的图像文件
%Z是原始图像z=f(x,y)的灰度值的矩阵
%Xnew,Ynew是在新的图像网格上对应矩阵的行列维数
%Znew是插值后新的图像矩阵
Z=f;%imread(f);
%imshow(Z);
%title('原始图像');
[dim(1) dim(2) dim(3)]=size(Z);%读取图像维数，用于判断是否彩色图像

Zold=double(Z);%转化为双精度型减小插值误差
for i=1:2
   for j=1:2
      Znew(i:2:2*dim(1),j:2:2*dim(2))=Zold;%上采样到二倍大小
   end
end
Xnew=2*dim(1);
Ynew=2*dim(2);
      
        for i=1:2:Xnew
            for j=2:2:Ynew
                P2=[Znew(i,max(j-3,1))    Znew(i,j-1)   Znew(i,min(j+1,Ynew-1))   Znew(i,min(j+3,Ynew-1))];
                a=[-(1/16) 9/16 9/16 -(1/16)]';
                Znew(i,j)=P2*a;
            end
        end
        for i=2:2:Xnew
            for j=1:2:Ynew
                P3=[Znew(max(i-3,1),j)    Znew(i-1,j)   Znew(min(i+1,Xnew-1),j)   Znew(min(i+3,Xnew-1),j)];
                a=[-(1/16) 9/16 9/16 -(1/16)]';
                Znew(i,j)=P3*a;
            end
        end
        for i=2:2:Xnew-2
            for j=2:2:Ynew-2
                P1 = [Znew(i,max(j-3,1))    Znew(i,j-1)   Znew(i,min(j+1,Ynew-1))   Znew(i,min(j+3,Ynew-1))];
                a=[-(1/16) 9/16 9/16 -(1/16)]';
                Znew(i,j)=P1*a;
               
            end
        end
        %后处理
for i=1:Xnew
    Znew(i,Ynew)=Znew(i,Ynew-1);
end
for j=1:Ynew
    Znew(Xnew,j)=Znew(Xnew-1,j);
end
Znew=uint8(Znew);
%figure;
%imshow(Znew);
%title('bicubic插值后图像');