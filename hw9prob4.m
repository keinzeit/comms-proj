%% EE107 HW9 Problem 4
% Kenny Yau

%% Part A
allData = [0 0 0 0
        0 0 0 1
        0 0 1 0
        0 0 1 1
        0 1 0 0
        0 1 0 1
        0 1 1 0
        0 1 1 1
        1 0 0 0
        1 0 0 1
        1 0 1 0
        1 0 1 1
        1 1 0 0
        1 1 0 1
        1 1 1 0
        1 1 1 1];
    
 allData(allData == 0) = -1 % should have made all the 0s -1s in the first place
 index = 0:15; % index for graphs

% Channel 1 (ak = 1)
ak = 1;
x1 = [0.1 -0.25 1 -0.25 0.1];

for i = 1:16
    data = allData(i,:);
    y1ak1(i) = channel(x1,data,ak);
end

figure(1)
subplot(1,2,1)
stem(index,y1ak1)
title('Channel 1 (ak = 1)')
xlabel('Decimal Representation of the 4-bit Data Sequence')
ylabel('Output Value')

% Channel 2
ak = 1;
x2 = [-0.2 0.5 1 0.5 -0.2];

for i = 1:16
    data = allData(i,:);
    y2ak1(i) = channel(x2,data,ak);
end

subplot(1,2,2)
stem(index,y2ak1)
title('Channel 2 (ak = 1)')
xlabel('Decimal Representation of the 4-bit Data Sequence')
ylabel('Output Value')

% Channel 1 (ak = -1
ak = -1;
x1 = [0.1 -0.25 1 -0.25 0.1];

for i = 1:16
    data = allData(i,:);
    y1akn1(i) = channel(x1,data,ak);
end

figure(2)
subplot(1,2,1)
stem(index,y1akn1)
title('Channel 1 (ak = -1)')
xlabel('Decimal Representation of the 4-bit Data Sequence')
ylabel('Output Value')

% Channel 2
ak = -1;
x2 = [-0.2 0.5 1 0.5 -0.2];

for i = 1:16
    data = allData(i,:);
    y2akn1(i) = channel(x2,data,ak);
end

subplot(1,2,2)
stem(index,y2akn1)
title('Channel 2 (ak = -1)')
xlabel('Decimal Representation of the 4-bit Data Sequence')
ylabel('Output Value')

%%
% From the graphs, we see that the output of channel 1 is always the same
% sign as the given ak, while for channel 2 there is a possibility that the
% sign of the output does not match that of ak. Thus, channel 2 may have
% detector errors due to ISI.

%% Part B
var = 0.1;
sigma = sqrt(var);
awgnoise = sigma*rand(1,16);
y1ak1noise = y1ak1 + awgnoise;
y2ak1noise = y2ak1 + awgnoise;
y1akn1noise = y1akn1 + awgnoise;
y2akn1noise = y2akn1 + awgnoise;


figure(3)
subplot(1,2,1)
stem(index,y1ak1noise)
title('Channel 2 (ak = 1) w/ noise')
xlabel('Decimal Representation of the 4-bit Data Sequence')
ylabel('Output Value')

subplot(1,2,2)
stem(index,y2ak1noise)
title('Channel 2 (ak = 1) w/ noise')
xlabel('Decimal Representation of the 4-bit Data Sequence')
ylabel('Output Value')

figure(4)
subplot(1,2,1)
stem(index,y1akn1noise)
title('Channel 2 (ak = -1) w/ noise')
xlabel('Decimal Representation of the 4-bit Data Sequence')
ylabel('Output Value')

subplot(1,2,2)
stem(index,y2akn1noise)
title('Channel 2 (ak = -1) w/ noise')
xlabel('Decimal Representation of the 4-bit Data Sequence')
ylabel('Output Value')

%%
% Adding white gaussian noise increases the probability of detector errors
% in channel 2 (we see more sign mismatches). It also affects channel 1, as
% the output values are closer to the threshold. Changing the variance to a
% higher value, I was able to get some detector errors in channel 1.

%% Part C
input = rand(1,128);
input(input>0.5) = 1;
input(input<=0.5) = -1;

% 32-sample pulse
T = 32;
t = linspace(0,T,32);
g = sin(pi/T*t);

%
output1 = conv(input, x1);
output2 = conv(input, x2);

figure(5)
for i = 1:length(output1)
    pulse1 = output1(i)*g;
    pulse2 = output2(i)*g;
    plot(pulse1,'b'); hold on
    plot(pulse2,'r'); hold on
end
title('Eye Diagram w/o Noise')
legend('Channel 1','Channel 2')

figure(6)
n2 = length(output2);
for i = 1:n2
    pulse1 = output1(i)*g + sigma*rand(1,32);
    pulse2 = output2(i)*g + sigma*rand(1,32);
    plot(pulse1,'b'); hold on
    plot(pulse2,'r'); hold on
end
title('Eye Diagram w/ Noise')
legend('Channel 1','Channel 2')

%%
% In the eye diagram w/o noise, the separation between pulse levels are
% pretty distinct, which means it is easy to discern the value of yk. When
% noise is added, that separation becomes really small, and values of one
% pulse cross over into the values that lower pulses can take. For a
% digital signal, this crossover is only a problem when the pulse amplitude
% is near the threshold, and we do see signs of the pulses with a positive
% no-noise value cross into the negative region, and vice versa.
