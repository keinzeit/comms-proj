%% Perfect Channel

% This sections leaves out all the modulation, channel, matched filter,
% equalization, and sampling and detection portions. Essentially it assumed
% a perfect channel - perfect transmission of the data.

filename = 'file.jpeg';
qbits = 16;

% Image Pre-processing
% r2*c2 represents the number of 8x8 blocks that fit in the image
[Zq2,r2,c2,m2,n2,minval2,maxval2] = ImagePreProcess_gray(filename,qbits);
% Conversion to Bit-stream
% R2 and C2 are the dimensions of the block - possibly unnecessary
[bs2,R2,C2,N2] = blocks2bitStream(qbits,Zq2);

% Conversion to Image
newZq2 = bitstream2blocks(bs2,qbits,R2,C2);
% Image Post-processing
ImagePostProcess_gray(newZq2,r2,c2,m2,n2,minval2,maxval2);

%%

% take 1-norm of the difference.0 means the newZq2 was constructed correctly
sum(sum(sum(abs(Zq2 - newZq2))))