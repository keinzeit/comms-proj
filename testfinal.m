%% final test
I0 = imread('file.jpeg');               % Import Image
I0 = imrotate(I0, 270);                 % Rotate Image
[Zt dctZ I] = ImPreProc(I0);            % Send image for processing

figure(1); imshow(I)                    % Show Processed Image        
figure(2); imshow(dctZ)                 % Show DCT of Image

qBits = 16;                             % Quantization bits
Zq = quantizer(qBits, Zt);              % Send DCT image to be quantized


% Create Bit Stream out of Quantized Data
[bitStream M N B] = Convert_to_Bitstream(qBits, Zq);

T = 32;                                 % Sampling Period
PAM_level = 2;                          % PAM levels

% Modulate Bit Stream using Half Sine Pulses
[modBsSP sPulse Am] = ModulatedSinePulse(bitStream, T, qBits, PAM_level, M, N, B);

figure(3); plot(sPulse)
figure(4); freqz(sPulse)

alpha = 0.5;                            % Roll- Off Factor
K = 4;

% Modulate Bit Stream using SRRC Pulses
[modBsSRRC modStreamSRRC] = ModulatedSRRCPulse(bitStream, alpha, T, K, Am, M, N, B, qBits, PAM_level);

% pulse time and frequency plots
figure(5); plot(modBsSRRC)
figure(6); freqz(modBsSRRC)
figure(7); plot(modBsSP)
figure(8); freqz(modBsSP)
figure(9); plot(modStreamSRRC)


