Author:

Xiangjun Zhang

Dept. of ECE, McMaster University
Hamilton, ON, Canada, L8R 3H5


Usage:

ARInterpolation.exe LR.pgm HR.pgm


This version can only support 8-bit gray scale images in pgm format. The zoom in factor is 2. The work is based on the paper:

Xiangjun Zhang and Xiaolin Wu, "Image Interpolation by Adaptive 2D Autoregressive Modeling and Soft-Decision Estimation", IEEE Trans. Image Processing, Vol. 17, No. 6, Page(s): 887-896, June 2008.


For other type of images, one can convert it to pgm file by using the included matlab code "convertPGM.m". An example of usage is: convertPGM('image1.tif', 'image1.pgm');


For further information, please contact: zhangxj@grads.ece.mcmaster.ca