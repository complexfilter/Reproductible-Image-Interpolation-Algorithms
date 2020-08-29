function PSNR = get_PSNR(image1, image2)
% Calculate the PSNR between two images of the same size. 
% If two images are identical, the PSNR is set to 0.
%
% Input:
% image1, image2: matrices of two images to compare. They should be of the
% same size.
%
% Output:
% PSNR: the PSNR between the two images.

% Guoshen Yu 23/02/05

image1 = double(image1);
image2 = double(image2);

[image1_length, image1_width] = size(image1);
[image2_length, image2_width] = size(image2);

if (image1_length ~= image2_length) || (image1_width ~= image2_width)
    error('Images should be of same size.');
end

dif = (image1 - image2).^2;
dif = mean(dif(:));

if dif == 0
    PSNR = 0;
else
    PSNR = 10 * log10((255^2)/dif);
end
        