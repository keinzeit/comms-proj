function [Received_Image] = Image_Post_Proc(bitStream, M, N, B, minZ, maxZ)

%% Function takes the bitstream obtained through sampling and converts it
%back to an image.

%Inputs bitStream
%Outputs image

%Convert binary stream back to decimal values

decimalBitStream = bi2de(bitStream);

newZ = reshape(decimalBitStream,[M,N,B]);

newZ = newZ./65536

dctZ = im2double(newZ*(maxZ-minZ)+minZ);

dctZ = reshape(dctZ, [M*log2(B) N*log2(B)]);

fun=@idct2;
newZ=blkproc(dctZ,[8 8],fun);
whos dctZ

figure()
imshow(newZ)






%fun=@idct2;
%iDCTZ=blkproc(scaledZ,[8 8],fun)



