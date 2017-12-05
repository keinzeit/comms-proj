function [bitStream M N B] = Convert_to_Bitstream(qBits, Zq)

%Function takes the 3d array of quantized matrices and converts them into a
%binary single stream of bits.  The row length of this bit stream is
%dependent on the number of quantization levels chosen

% inputs Zq
% outputs bitStream

[M,N,B] = size(Zq);
 
len = M*N;
bsd = zeros(M*N*B,1);                       %preallocate memory for vector
for i = 1:B
    bs_temp = reshape(Zq(:,:,i),len,1);     % reshape B MxN dct blocks into a long column vector
    start = len*(i-1)+1;
    stop  = len*i;
    bsd(start:stop) = bs_temp;              % put values into correct location in bs
end
 
bitStream = de2bi(bsd,qBits);               % convert quantized values into q-bit numbers. right MSB; use correct flag for left MS

