function y=sri1(x)

[M,N]=size(x);
%e=zeros(M,N);
for i=1:2
   for j=1:2
      y(i:2:2*M,j:2:2*N)=x;
   end
end

T=5;cnt=1;
ix=[-1 -1 1 1];
iy=[-1 1 -1 1];
mx=[0 0 1 1];
my=[0 1 0 1];
cnt=1;
for i=-2:4
   for j=-2:4
      nx(cnt)=i;ny(cnt)=j;cnt=cnt+1;
  end
end
%nx=[0 1 1 0 -1 -1  0  1 2 2 0 1 -1 -1 2 2];
%ny=[0 0 1 1  0  1 -1 -1 0 1 2 2 -1  2 -1 2];

th=8;
for i=T:M-T
   for j=T:N-T
      for k=1:4
      C(:,k)=diag(x(i+ix(k)+nx,j+iy(k)+ny));
   	end
      r=diag(x(i+nx,j+ny));
      s=diag(x(i+mx,j+my));
      if det(C'*C)==0|var(s)<th
         a=ones(4,1)/4;
      else
         a=inv(C'*C)*(C'*r);
         %a=fun(C,r);
      end
      y(2*i,2*j)=sum(a.*s); 
      if y(2*i,2*j)<0|y(2*i,2*j)>255
         a=ones(4,1)/4;y(2*i,2*j)=sum(a.*s);
		end
   end
end


