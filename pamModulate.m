function modSignal = pamModulate(bitStream,pulse,spb,PAM_level)
% function modSignal = pamModulate(bitStream,pulse,spb,PAM_level)
% Function takes in a bitStream and performs PAM Modulation with the
% supplied pulse shaping filter.

if nargin < 4
    PAM_level = 2;
end

R = log2(PAM_level);    % number of bits used to change amplitude of pulse shape
                        % so far, unused. in future, code may be rewritten
                        % so that the function can perform more than just
                        % binary PAM
                        
Am(bitStream==1) = 1;
Am(bitStream==0) = -1;


L = length(pulse);

modSignal = zeros(1,(length(Am)+L-1)*spb+1); 
tModSignal = (0:length(modSignal)-1)/spb; % not needed for now
for i = 1:length(Am)
    start = (i-1)*spb+1;
    stop  = (i-1)*spb+L;
    modSignal(start:stop) = modSignal(start:stop) + Am(i)*pulse; %consider removing modStreamHS for half-sine pulse
end

return
% % SRRC Pulse specific
% modStreamSRRC = zeros(1,(size(Am,1)*size(Am,2)+2*K-1)*T+1);
% % tModBsSP = (1:length(modStreamSRRC))/T;
% for i = 1:length(Am)
%     start = (i-1)*spb+1;
%     stop  = (i-1)*spb+L;
%     modStreamSRRC(start:stop) = modStreamSRRC(start:stop) + Am*SRRCPulse;
% end
