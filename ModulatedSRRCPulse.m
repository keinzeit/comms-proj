function [SRRCPulse,modStreamSRRC] = ModulatedSRRCPulse(bitStream,alpha,T,K,Am,M,N,B,qBits,PAM_level)
% Function inputs the bitStream and modulates the bit stream using a Square
% Root Raised Cosine shaping filter.  The function creates the SRRC pulse
% and then modulates the bit stream

t2 = linspace(-K*T,K*T,2*K*T+1);
R = log2(PAM_level); % number of bits used to change amplitude of pulse shape
pulsesPerRow = qBits/R; % number of pulses made from a row of Am

%Create the SRRC pulse
SRRCPulse = zeros(1, length(t2));
for i = 1:length(SRRCPulse)
    if t2(i) == 0
        %denominator equals 0
        SRRCPulse(i) = 1 - alpha + 4*alpha/pi;
    elseif abs(t2(i)) == T/(4*alpha)
        %denominator equals 0
        SRRCPulse(i) = (alpha/sqrt(2)) * ((1+(2/pi))*sin(pi/(4*alpha)) + (1-(2/pi))*cos(pi/(4*alpha)));
    else
        SRRCPulse(i) = (sin(pi*t2(i)/T*(1-alpha)) + 4*alpha*t2(i)/T*cos(pi*t2(i)/T*(1+alpha))) / (pi*t2(i)/T*(1-((4*alpha*t2(i)/T)^2)));
    end 
end


modStreamSRRC = zeros(1,(size(Am,1)*size(Am,2)+2*K-1)*T+1);
% tModBsSP = (1:length(modStreamSRRC))/T;
for i = 1:M*N*B*qBits
    j = fix((i-1)/pulsesPerRow)+1;
    k = mod(i-1,pulsesPerRow)+1;
    start = (i-1)*T+1;
    stop  = (i-1)*T+2*K*T+1;
    modStreamSRRC(start:stop) = modStreamSRRC(start:stop) + Am(j,k)*SRRCPulse; %consider removing modBsSP for half-sine pulse  
end


