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

figure(113)
plot(tModBsSP(1:T*10),modBsSP(1:T*10))
title('10 Random Symbols w/ Half-Sine Pulse')
xlabel('Time (s)')
ylabel('Amplitude')

figure(114)
plot(

freqz(tModBsSP)

freqz(mod

delay = zeros(1, 31);                           %vector of 31 zeros to space between values in h
h = [1 delay 1/2 delay 3/4 delay -2/7];         % channel filter FIR

figure(115)
freqz(h)

figure(116)
stem(h)
axis([-10 105 -.5 1.5])

ysine = conv(modBsSP(1:T*10), h);                %channel output 
figure(117)
plot(ysine)


ysrrc = conv(modStreamSRRC(1:T*10), h);%channel output

eyediagram(ysine,32)

eyediagram(ysrrc,32)

figure(118)
plot(ysrrc)


figure(119)
plot(ysn)

figure(120)
eyediagram(ysn,32)

figure(121)
eyediagram(ysrrcn,32)
