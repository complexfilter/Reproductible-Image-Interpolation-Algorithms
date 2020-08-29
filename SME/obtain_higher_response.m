load qmf.mat;
qmf_filter=qmf;
tttt=Wavelet_lf(:,:,1);
llll=cconv2(qmf_filter(end:-1:1),qmf_filter(end:-1:1),imL,-5,-4);
figure;imagesc(llll-tttt);

llll=cconv2(qmf_filter,qmf_filter,imL,-5,-5);
figure;imagesc(llll-tttt);


ttt=qmf_filter(end:-1:1);

%ggg1=firpm(1350,[0 0.495 0.505 1],[1 1 0 0],[200 1]);
ggg1=firpm(2100,[0 0.500 0.505 1],[1 1 0 0],[200 1]);
ggg2=firpm(510,[0 0.5 0.53 1],[1 1 0 0],[200 1]);
ggg2=firpm(2099,[0 0.500 0.505 1],[1 1 0 0],[200 1]);
fvtool(ggg1,1,ggg2,1,upsample(ttt,2,0),1);
final=2*conv(conv(conv(ggg1,ggg1),ggg1),upsample(ttt,2,0));
fvtool(final(3160-2*5-1:2:3160+5*2),1,ttt,1);
fvtool(final(3160-2*5:2:3160+5*2),1,ttt,1);
fvtool(final(3160-4*5:1:3160+5*3),1,ttt,1);
% 

ggg1=firpm(1350,[0 0.495 0.505 1],[1 1 0 0],[200 1]);
ggg2=firpm(510,[0 0.5 0.53 1],[1 1 0 0],[200 1]);
final=2*conv(conv(ggg1,ggg2),upsample(ttt,2,0));

fvtool(final(941-2*5-16:941+2*2+16))
% 
qmf_2x=final(941-2*5-16:941+2*2+16);
% 
% fvtool(qmf_2x(2:2:end))

ddd=1/4*cconv2(qmf_2x,qmf_2x,I_ini,-25,-23);

figure;imagesc(ddd(1:2:end-1,1:2:end-1)-llll);
figure;imagesc(ddd(2:2:end-1,2:2:end-1)-llll);

% 
% 
% fvtool(final(3155-10*5:2:3155+10*2),1,ttt,1);
% 
% l=fft(ttt);
% 
% dd=ifft([l(1:4) 0 0 0 0 0 0 0 0 l(5:8)]);

% I am going to do it anyway.

ll=cconv2(upsample(ttt,2,0),upsample(ttt,2,0),cconv2(upsample(ttt,2,0),upsample(ttt,2,0),I_ini,-9,-9)/2,-9,-9)/2;
figure;imagesc(ll(3:end,3:end)-imH_linear_lf);


ll=cconv2(qmf_2x,upsample(ttt,2,0),I_ini,-9,-9)/2;
figure;imagesc(ll(3:end,3:end)-imH_linear_lf);

ggg=2*upsample(ttt,2,0);
ggg=[zeros(1,16) ggg zeros(1,15)];
qmf_2x_higher=ggg-qmf_2x;

ddd1=1/4*cconv2(qmf_2x,qmf_2x,f,-25,-23);
ddd2=1/4*cconv2(qmf_2x,qmf_2x_higher,f,-25,-23);
ddd3=1/4*cconv2(qmf_2x_higher,qmf_2x,f,-25,-23);
ddd4=1/4*cconv2(qmf_2x_higher,qmf_2x_higher,f,-25,-23);
ddd=ddd1+ddd2+ddd3+ddd4;
figure;imagesc(ddd(2:2:end,2:2:end)-llll);


