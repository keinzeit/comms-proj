function [modStreamHS,sPulse,Am] = ModulatedSinePulse(bitStream,T,qBits,PAM_level,M,N,B)
% Function inputs the bit stream, period, and PAM level.  Using this data, 
% the function will create a Half-Sine pulse and modulate the bit stream with the pulse
% to prepare thee data for the channel

%inputs bitStream, T, qBits, PAM_level
%outputs Half_Sine_bsMod

% Creation of Half-Sine pulse

t1 = linspace(0,T,T+1);                 % Pulses include the first value of the next pulse
sPulse = sin(pi*t1/T);

R = log2(PAM_level);                    % number of bits used to change amplitude of pulse shape
pulsesPerRow = qBits/R;                 % number of pulses made from a row of Am

Am = zeros(M*N*B,R);
for i = 1:M*N*B
    for j = 1:pulsesPerRow 
        head = 1+log2(PAM_level)*(j-1);
        tail = log2(PAM_level)*j;
        Am(i,j) = bi2de(bitStream(i,head:tail));
    end
end
if PAM_level == 2
    Am(Am==0) = -1;
end

modStreamHS = zeros(1,size(Am,1)*size(Am,2)*T+1); 
tModBsSP = (1:length(modStreamHS))/T;
for i = 1:M*N*B*qBits
    j = fix((i-1)/pulsesPerRow)+1;
    k = mod(i-1,pulsesPerRow)+1;
    start = (i-1)*T+1;
    stop  = i*T+1;
    modStreamHS(start:stop) = Am(j,k)*sPulse; %consider removing modStreamHS for half-sine pulse
    
end


