function [Zq] = quantizer(qBits, Zt)

% Function quantizes the 8x8xB matrix into 2^qbits levels.  quantizer is
% pre-set to accept qbits = 8,16

% Inputs Zt, qBits
% Outputs Zq

if qBits == 8
   Zq=im2uint8(Zt);                    % Quantize to 2^8 levels
elseif qBits == 16
   Zq=im2uint16(Zt);                   % Quantize to 2^16 levels
end       

whos Zq