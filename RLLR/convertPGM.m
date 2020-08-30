function convertPGM(inputfile, outputfile)

img = imread(inputfile);
if max(max(img))<256
    maxval = 255;
    flag = 0;
else
    maxval = 65535;
    flag = 1;
end

fid = fopen(outputfile, 'wb');
fprintf(fid, 'P5\n%d %d\n%d\n', size(img,2), size(img,1), maxval);
if flag==0
    fwrite(fid, img', 'uint8');
else
    fwrite(fid, img', 'uint16');
end
fclose(fid);