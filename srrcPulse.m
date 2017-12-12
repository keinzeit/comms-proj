function SRRCPulse = srrcPulse(alpha,spb,K)
% function SRRCPulse = srrcPulse(alpha,spb,K,PAM_level)
% Function returns the impulse response of a square-root raised cosine
% filter. 

t = linspace(-K*spb,K*spb-1,2*K*spb);

%Create the SRRC pulse
SRRCPulse = zeros(1, length(t));
for i = 1:length(SRRCPulse)
    if t(i) == 0
        %denominator equals 0
        SRRCPulse(i) = 1 - alpha + 4*alpha/pi;
    elseif abs(t(i)) == spb/(4*alpha)
        %denominator equals 0
        SRRCPulse(i) = (alpha/sqrt(2)) * ((1+(2/pi))*sin(pi/(4*alpha)) + (1-(2/pi))*cos(pi/(4*alpha)));
    else
        SRRCPulse(i) = (sin(pi*t(i)/spb*(1-alpha)) + 4*alpha*t(i)/spb*cos(pi*t(i)/spb*(1+alpha))) / (pi*t(i)/spb*(1-((4*alpha*t(i)/spb)^2)));
    end 
end

return