%% EE107 Final Project
% Kenny Yau, Jon Manninen
%
% This project is a simulation of the model of communication systems for
% images, such as that used for digital HDTV. It consists of a number of
% parts:
%     - Image pre-processing
%     - Conversion to bit-stream
%     - Modulation
%     - Channel Effects
%     - Noise
%     - Matched Filter
%     - Receiver Equalization
%     - Sampling and Detection
%     - Converion to Image
%     - Image post-processing

% Stuff to do:
Image pre-process:
Find own image
convert image to double data type using im2double
make sure image is m-by-n, where both m and n are divisible by 8
perform DCT image in 8x8 blocks using blkproc
scale all values to fit between 0 and 1
quantize the scaled, transformed image using code from hw 7
use reshape and permute to convert the scaled image? into what?
    - reshape reorganizes a matrix with certain dimensions and size into a matrix of another dimension and size
    - we are organizing the image data into 8x8 blocks by stacking them on top of each other
    - use permute to make sure first two dimensions specify a single 8x8 block

Conversion to bit-stream:
we transmit N DCT blocks at a time
Convert N DCT blocks into a single column vector
convert each element in this column vector into an 8-bit binary number using de2bi, stored as a row
    - ex [5] -> [1 0 1 0 0 0 0 0]

Modulation:
multiply bit value by 32 samples of half-sine pulse and 2K*32 samples of SRRC pulse
 - SRRC pulse lasts longer than 1 bit duration, there will be overlap of 2K-1 values
 - normalize SRRC by a factor A such that power of both pulses are the same
need to plot eye diagram and modulated signal!! plus both pulse-shaping functions

Channel:
represent channel h[n] as vector, make sure it is power of 2, make sure each delta occurs at start of next bit
send modulated data through channel (use conv on each row of the modulated signal matrix)
plot freq response and eye diagram

Noise:
add noise to output of channel
plot eye diagrams of both pulse-shaping functions after noise

Matched Filter:
come up with the matched filter for each pulse shape
plot impulse and frequency response of both MFs
plot eye diagram at output of each MF

Equalizer:
Zero-Forcing (ZF): is the inverse of the channel response
use filter to make ZF
MMSE


