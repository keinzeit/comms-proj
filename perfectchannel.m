%% Perfect Channel

% This sections leaves out all the modulation, channel, matched filter,
% equalization, and sampling and detection portions. Essentially it assumed
% a perfect channel - perfect transmission of the data.

filename = 'file.jpeg';
qbits = 16;

% Image Pre-processing
% r2*c2 represents the number of 8x8 blocks that fit in the image
[Zq2,r2,c2,m2,n2,minval2,maxval2] = ImagePreProcess_gray(filename,qbits);

% Conversion to Bit-stream
% R2 and C2 are the dimensions of the block - possibly unnecessary
[bs2,R2,C2,N2] = blocks2bitStream(qbits,Zq2);

% % Modulation w/ HS Pulse
% PAM_level = 2; % PAM levels
% spb = 32;      % Samples per bit
% t1 = linspace(0,spb-1,spb); % Pulse does not include first value in the next pulse
% HSPulse = sin(pi*t1/spb);
% % Modulate Bit Stream using Half Sine Pulses
% HS_mod = pamModulate(bs2,HSPulse,spb);
% figure
% plot(HS_mod(1:1000))
% 
% % Modulation w/ SRRC Pulse
% alpha = 0.5;
% K = 2;
% SRRCPulse = srrcPulse(alpha,spb,K);
% figure
% plot(SRRCPulse)
% SRRC_mod = pamModulate(bs2,SRRCPulse,spb);
% figure
% plot(SRRC_mod(1:2^10))

% Conversion to Image
newZq2 = bitstream2blocks(bs2,qbits,R2,C2);
% Image Post-processing
newZ = ImagePostProcess_gray(newZq2,r2,c2,m2,n2,minval2,maxval2);
figure
imshow(newZ)
%%
% take 1-norm of the difference.0 means the newZq2 was constructed correctly
sum(sum(sum(abs(Zq2 - newZq2))));