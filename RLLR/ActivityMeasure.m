
[height,width,channel] = size(image);
y = image(:,:,channel);
image_size = height*width;
detail_mask = zeros(image_size);

MP = 14;
Eim = [repmat(y(1,:),MP,1);y;repmat(y(m,:),MP,1)];
Eim = [repmat(Eim(:,1),1,MP),Eim,repmat(Eim(:,n),1,MP)];
[M N]=size(Eim);%MÐÐNÁÐ
my=[-1 -1  1 1];
mx=[-1  1 -1 1];

for i=1+MP:M-MP
   for j=1+MP:N-MP
    neighbor = diag(Eim(i+my,j+mx));
    if var(neighbor)<th
       masked_image((i-MP),(j-MP))=255;  
       count = count +1 ;
    end
   end
end

function detail_mask = ActivityMeasure(neighbor)
neighbor_size = size(neighbor);
histogram = zeros(256);
for i  = 1:neighbor_size
    histogram(neighbor(i)) = histogram(neighbor(i))+1;
end
