load qmf_2x.mat;
load qmf.mat
qmf_filter=qmf;
ttt=qmf(end:-1:1);
ggg=2*upsample(ttt,2,0);
ggg=[zeros(1,16) ggg zeros(1,15)];
qmf_2x_higher=ggg-qmf_2x;

ddd1=1/4*cconv2(qmf_2x,qmf_2x,f,-27,-27);
ddd2=1/4*cconv2(qmf_2x,qmf_2x_higher,f,-27,-27);
ddd3=1/4*cconv2(qmf_2x_higher,qmf_2x,f,-27,-27);
ddd4=1/4*cconv2(qmf_2x_higher,qmf_2x_higher,f,-27,-27);
ddd=ddd1+ddd2+ddd3+ddd4;
llll=cconv2(qmf_filter(end:-1:1),qmf_filter(end:-1:1),imL,-6,-6);
%% There exists some difference between llll and  Wavelet(:,:,1).
distance=llll-Wavelet(:,:,1);
figure;imagesc(abs(distance(50:200,50:200)));

figure;imagesc(abs(llll-Wavelet(:,:,1)));
 

figure;imagesc(ddd(2:2:end,2:2:end)-llll);
figure;imagesc(abs(llll));

ddd1_dsp2=ddd1(2:2:end,2:2:end);
ddd2_dsp2=ddd2(2:2:end,2:2:end);
ddd3_dsp2=ddd3(2:2:end,2:2:end);
ddd4_dsp2=ddd4(2:2:end,2:2:end);

lll1=1/4*cconv2(qmf_2x       ,qmf_2x       ,upsample(upsample(ddd1_dsp2,2,1).',2,1).',-25,-23);
lll2=1/4*cconv2(qmf_2x       ,qmf_2x_higher,upsample(upsample(ddd2_dsp2,2,1).',2,1).',-25,-23);
lll3=1/4*cconv2(qmf_2x_higher,qmf_2x       ,upsample(upsample(ddd3_dsp2,2,1).',2,1).',-25,-23);
lll4=1/4*cconv2(qmf_2x_higher,qmf_2x_higher,upsample(upsample(ddd4_dsp2,2,1).',2,1).',-25,-23);

figure;imshow(uint8(lll1+lll2+lll3+lll4));

norm(ddd1_dsp2(50:200,50:200)-llll(50:200,50:200),'fro')/norm(llll(50:200,50:200),'fro')

figure;imagesc((ddd1_dsp2(50:200,50:200)));
figure;imagesc((ddd2_dsp2(50:200,50:200)));
figure;imagesc((ddd3_dsp2(50:200,50:200)));
figure;imagesc((ddd4_dsp2(50:200,50:200)));