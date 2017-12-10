function Qmmse = qmmse(channel,N,noise)
% function Qmmse = qmmse(channel, noise)
% Returns the MMSE Equalizer for a given channel. Noise power is optional.

% set default values
if nargin < 2
    N = 65536; % this value was chosen according to the length of our output vector
end
if nargin < 3
    noise = 0;
end

channel_FFT = fft(channel,N);

% create MMSE Equalizer
Qmmse = (conj(channel_FFT)./((abs(channel_FFT)).^2+2*noise));

return
