function [bitStream,R,C,N] = blocks2bitStream(qBits,Zq)
% Function takes the 3d array of quantized matrices and converts them into a
% binary single stream of bits.  The row length of this bit stream is
% dependent on the number of quantization levels chosen

% inputs Zq
% outputs bitStream

[R,C,N] = size(Zq);
 
len = R*C;
bsd = zeros(R*C*N,1);                       %preallocate memory for vector
for i = 1:N
    bs_temp = reshape(Zq(:,:,i),len,1);     % reshape N RxC dct blocks into a long column vector
    start = len*(i-1)+1;
    stop  = len*i;
    bsd(start:stop) = bs_temp;              % put values into correct location in bs
end

bitArray = de2bi(bsd,qBits);               % convert quantized values into q-bit numbers. right MSB; use correct flag for left MS

bitStream = array2stream(bitArray);

% % used to test my Rx side code
% RxBitArray = stream2array(bitStream,qBits);
% diffArray = bitArray - RxBitArray;
% sum(sum(abs(diffArray))) % 0 means received array is in same order as original

return

