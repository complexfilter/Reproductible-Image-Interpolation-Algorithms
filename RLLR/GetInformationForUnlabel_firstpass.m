function [S,Y_hat] = GetInformationForUnlabel_firstpass(i,j,win,swin,Eim,options)
     [M,N] = size(Eim);
      k=1;
      kernel = options.kernel;
      for ii=i-2*(win-1):2:i-1
          for jj=j-2*(win-1):2:j+2*win-1
              pixel_no = ii*N+jj;
              sample = Eim(ii-2*swin:2:ii+2*swin,jj-2*swin:2:jj+2*swin);
              sample_kernel = sample.*kernel;
              sample_kernel = sample_kernel(:);
              X(:,k) = sample_kernel;                   
              Y_hat(1,k) = Eim(ii,jj);
              k=k+1;              
          end
      end      
         
      ii = i;
      for jj=j-2*(win-1):2:j-2
          pixel_no = ii*N+jj;
          sample = Eim(ii-2*swin:2:ii+2*swin,jj-2*swin:2:jj+2*swin);
          sample_kernel = sample.*kernel;
          sample_kernel = sample_kernel(:);
          X(:,k) = sample_kernel;
          Y_hat(1,k) = Eim(ii,jj);
          k=k+1; 
      end
      
      center = Eim(i-swin:i+swin,j-swin:j+swin);
      CB_kernel = center.*options.kernel;
      CB_kernel = CB_kernel(:);
      CB_kernel = repmat(CB_kernel,1,k-1);
      
      S=Computweight_matrix(CB_kernel,X);

end
