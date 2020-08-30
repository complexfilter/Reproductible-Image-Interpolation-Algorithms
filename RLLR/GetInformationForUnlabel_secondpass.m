function [S,Y_hat] = GetInformationForUnlabel_secondpass(i,j,Xl,win,swin,Eim,flag,options)
      Xu_x = reshape(Xl(:,2),2*win,2*win);
      Xu_sub_x = Xu_x(1:7,1:7);
      Xu_x = reshape(Xu_sub_x,49,1);
      Xu_y = reshape(Xl(:,1),2*win,2*win);
      Xu_y = Xu_y+1;
      Xu_sub_y = Xu_y(1:7,1:7);
      Xu_y = reshape(Xu_sub_y,49,1);
      Xu(:,1) = Xu_y;
      Xu(:,2) = Xu_x;
      dim = size(Xu);
      [M,N] = size(Eim);
      k = 1;
      
      if flag == 0
          number = 23;
      else
          number = 25;
      end

      for i = 1:number
          ii = Xu(i,1);
          jj = Xu(i,2);
          pixel_no = ii*N+jj;
          sample = Eim(ii-2*swin:2:ii+2*swin,jj-2*swin:2:jj+2*swin);
          sample_kernel = sample.*options.kernel;
          sample_kernel = sample_kernel(:);
          X(:,k) = sample_kernel;
          Y_hat(1,k) = Eim(ii,jj);
          k=k+1; 
      end
      
      center = Eim(i-swin:i+swin,j-swin:j+swin);
      CB_kernel = center.*options.kernel;
      CB_kernel = CB_kernel(:);
      CB_kernel = repmat(CB_kernel,1,number);%÷ÿ∏¥L±È
      

      S = Computweight_matrix(CB_kernel,X);
   
end
