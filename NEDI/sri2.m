% second step of edge-directed interpolation
function y=sri2(x)
y=x;
[M,N]=size(x);
T=8;cnt=1;
mx=[0 0 1 -1];
my=[-1 1 0 0];
%nx=[-1 0 0 1 -2 1 -1 2 1 2 -2 -1 3 -3 0 0]; 
%ny=[0 -1 1 0 1 -2 2 -1 2 1 -1 -2 0 0 -3 3];
cnt=1;
for i=-2:3
   for j=-2:3
      nx(cnt)=i+j-1;ny(cnt)=i-j;cnt=cnt+1;
   end
end

th=8;
for i=T:2:M-T
   for j=T+1:2:N-T
       for k=1:4
           C(:,k)=diag(x(i+2*mx(k)+nx,j+2*my(k)+ny));
       end
      r=diag(x(i+nx,j+ny));
      s=diag(x(i+mx,j+my));
      if det(C'*C)==0|var(s)<th
         a=ones(4,1)/4;
     else
          a=inv(C'*C)*(C'*r);
%          %a=fun(C,r);
end
      y(i,j)=sum(a.*s);
      if y(i,j)<0|y(i,j)>255
         a=ones(4,1)/4;y(i,j)=sum(a.*s);
		end
   end
end

for i=T+1:2:M-T
   for j=T:2:N-T
       for k=1:4
       C(:,k)=diag(x(i+2*mx(k)+nx,j+2*my(k)+ny));
   end
       r=diag(x(i+nx,j+ny));
      s=diag(x(i+mx,j+my));
      if det(C'*C)==0|var(s)<th
         a=ones(4,1)/4;
     else
         a=inv(C'*C)*(C'*r);
%          %a=fun(C,r);
end
      y(i,j)=sum(a.*s);
      if y(i,j)<0|y(i,j)>255
         a=ones(4,1)/4;y(i,j)=sum(a.*s);
		end
   end
end


