function bitStream = array2stream(bitArray)
% function bitStream = array2stream(bitArray)

[R,C] = size(bitArray);
bitArray = bitArray';
bitStream = reshape(bitArray,[1,R*C]);

return