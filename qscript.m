% run this to answer the questions for the project
% assumes you have run finalproj.m to get the right variables in the
% workspace

% half-sine pulse
t1 = linspace(0,T,T+1); % I have pulses include the first value of the next pulse
t11 = linspace(0,1,T+1);
sPulse = sin(pi*t1/T); 

figure(101); plot(t11,sPulse)
title('half-sine pulse'); xlabel('Time (s)')
figure(102); freqz(sPulse)

% srrc pulse
A = norm(sPulse,2);
a = 0.5; % roll-off factor
t2 = linspace(-K*T,K*T,2*K*T+1);
t21 = linspace(-K,K,2*K*T+1);
srrcPulse = A*rcosdesign(a,K,2*T,'sqrt'); % energy of rcosdesign is one, just multiply by energy of half-sine pulse

figure(111); plot(t21,srrcPulse)
figure(112); freqz(srrcPulse)