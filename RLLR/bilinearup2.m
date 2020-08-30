function im2=bilinearup2(im)
im=double(im);
Xstep=0.5;
Ystep=0.5;
[m,n]=size(im);
for i=1:2
   for j=1:2
      im2(i:2:2*m,j:2:2*n)=im;%�ϲ�����������С
   end
end

for i=1:2:2*m
    for j=2:2:2*n-1
        im2(i,j)=0.5*(im2(i,j-1)+im2(i,j+1));
    end
end
for j=1:2:2*n
    for i=2:2:2*m-1
        im2(i,j)=0.5*(im2(i-1,j)+im2(i+1,j));
    end
end
a=ones(1,4)/4;
for i=2:2:2*m-1
    for j=2:2:2*n-1
        sam=[im2(i-1,j-1) im2(i-1,j+1) im2(i+1,j-1) im2(i+1,j+1)];
        im2(i,j)=sum(sam.*a);
    end
end
end