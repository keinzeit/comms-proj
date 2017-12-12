function blocks = bitstream2blocks(bitStream,qbits,R,C)
% function blocks = bitstream2blocks(bitStream)

bitArray = stream2array(bitStream,qbits);
decimals = bi2de(bitArray);

% convert double to uint
if qbits == 8
   decimals=uint8(decimals);
elseif qbits == 16
   decimals=uint16(decimals);
end       

len = length(decimals);
blocks = reshape(decimals,[R,C,len/(R*C)]);

return