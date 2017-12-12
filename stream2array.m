function bitArray = stream2array(bitStream,qbits)
% function bitArray = stream2array(bitStream)
% returns an array of bits, where each row of the array corresponds to 1
% number
% ex. [ 1 0 1 1  represents [11
%       0 1 0 1]              5]

numVals = length(bitStream)/qbits;
bitArray = reshape(bitStream,[qbits,numVals])';

return