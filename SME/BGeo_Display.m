function QQ_crop=BGeo_Display(BGeo,index)
QQ=zeros(256);
tt=BGeo(:,index);
QQ(tt(tt>0))=1;
QQ_crop=QQ(1:12,1:12);
end