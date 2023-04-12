% EE-442 Wireless Receivers: algorithms and architectures
% Final Project: OFDM Audio Transmission System
% EPFL Fall 2022
% Authors: Dong Chu, Han Miao, Huan-Ying Yeh

%% simple_ofdm.m
%
%% [Important] This code is modified from the original implementation from
% The outputs are only used for initial testing and sanity check. The project
% did not use this code and only discusses code that implemented ourselves.
%
%% Credit:
% This simulation code is based on "Jordan Street's" OFDM simulation
% presentation 'https://www.youtube.com/watch?v=SyKJrrNhPO8'. 

% Authors:
% 1- Naveed Ahmed: 'https://www.ntnu.edu/employees/naveed.ahmed'
% 2- Shujaat Khan:  Shujaat123@gmail.com

%% Implements a whole OFDM transmission system without continuous phase tracking (QPSK)
function [X, X_hat_none, X_hat, rec_im_none, rec_im, ber_none, ber] = simple_ofdm(im, snr, n_taps, ray_channel)

% Add AWGN and flat fading of desired SNR and #taps.
%
% Params:
% im - input image file
% training_type - 'block' or 'comb' or 'viterbi'
% snr - SNR
% n_taps - fading tap numbers
% ray_channel - [Optional] the ray_channel object
%
% Outputs:
% - X: transmitted bits from the image
% - X_hat_none: received bits
% - X_hat: received symbols with phast tracking [Result not used in report]
% - rec_im_none: the received image without phase tracking
% - rec_im: received image with block-type phase tracking [Result not used in report]
% - ber_none: BER without phase tracking
% - ber: BER with block-type phase tracking [Result not used in report]

%fft Size
nfft = 64;
n_fft = 64;

% Size of cycle prefix extension
n_cpe = 32;

% Number of channel taps.
mod_order = 2;

%% Read the image and convert it into binary format.
im_bin = dec2bin(im(:))';
im_bin = im_bin(:);

%% Binary stream to symbols
% 1. parse binary stream into mod_order bit symbols
% 2. pads input signal to appropriate length
sym_rem = mod(mod_order-mod(length(im_bin),mod_order),mod_order);
padding = repmat('0',sym_rem,1);
im_bin_padded = [im_bin;padding];
cons_data = reshape(im_bin_padded,mod_order,length(im_bin_padded)/mod_order)';
cons_sym_id = bin2dec(cons_data);

%% Symbol modulation
mod_ind = 2^(mod_order-1);
n = 0:pi/mod_ind:2*pi-pi/mod_ind;
in_phase = cos(n+pi/4);
quadrature = sin(n+pi/4);
symbol_book = (in_phase + quadrature*1i);
X = symbol_book(cons_sym_id+1);

%% IFFT move to time domain
% pad input signal to appropriate length
fft_rem = mod(n_fft-mod(length(X),n_fft),n_fft);
X_padded = [X, zeros(fft_rem,1)'];
X_blocks = reshape(X_padded,nfft,length(X_padded)/nfft);
x = ifft(X_blocks);

% Add cyclic prefix entension and shift from parallel to serial
x_cpe = [x(end-n_cpe+1:end,:);x];
x_s   = x_cpe(:);

%% Add AWGN
% Calculate data power
data_pwr = mean(abs(x_s.^2));

% Add noise to the channel
noise_pwr = data_pwr/10^(snr/10);
noise = normrnd(0,sqrt(noise_pwr/2),size(x_s))+normrnd(0,sqrt(noise_pwr/2),size(x_s))*1i;
x_s_noise = x_s + noise;

%% Apply fading channel 
g = exp(-(0:n_taps-1));
g = g/norm(g);
x_s_noise_fading = conv(x_s_noise,g,'same');

%% Deal with Rayleigh Channels if applicable
if strcmp(ray_channel, 'true') % Then a ray channel is supplied
    x_s_noise_fading = ray_channel(x_s_noise);
end

%% FFT move to frequency domain
% Remove cyclic prefix extension and shift from serial to parallel
x_p = reshape(x_s_noise_fading,nfft+n_cpe,length(x_s_noise_fading)/(nfft+n_cpe));
x_p_cpr = x_p(n_cpe+1:end,:);

% Move to frequency domain
X_hat_blocks = fft(x_p_cpr);

X_hat_none = X_hat_blocks(:);
X_hat_none = X_hat_none(1:end-fft_rem);

%% Estimate channels
G = X_hat_blocks(:,1)./X_blocks(:,1);
X_hat_blocks = X_hat_blocks./repmat(G,1,size(X_hat_blocks,2));

%% Symbol demodulation
% remove fft padding 
X_hat = X_hat_blocks(:);
X_hat = X_hat(1:end-fft_rem);

% Recover data from modulated symbols
A=[real(symbol_book) imag(symbol_book)];
if (size(A,2)>2)
    A=[real(symbol_book)' imag(symbol_book)'];
end

% Recover symbols for with or without estimations
rec_syms_est = knnsearch(A,[real(X_hat) imag(X_hat)])-1;
rec_syms_none = knnsearch(A,[real(X_hat_none) imag(X_hat_none)])-1;

% Parse to binary stream to remove symbol padding
rec_syms_cons = dec2bin(rec_syms_est);
rec_im_bin = reshape(rec_syms_cons',numel(rec_syms_cons),1);
rec_im_bin = rec_im_bin(1:end-sym_rem);
ber = sum(abs(rec_im_bin-im_bin))/length(im_bin);

%% Recover image (with channel estimation)
% rec_im = reshape(rec_im_bin,9,numel(rec_im_bin)/8);
rec_im = reshape(rec_im_bin,8,numel(rec_im_bin)/8);
rec_im = uint8(bin2dec(rec_im'));
rec_im = reshape(rec_im,size(im));

% Parse to binary stream to remove symbol padding
rec_syms_cons_none = dec2bin(rec_syms_none);
rec_im_bin_none = reshape(rec_syms_cons_none',numel(rec_syms_cons_none),1);
rec_im_bin_none = rec_im_bin_none(1:end-sym_rem);
ber_none = sum(abs(rec_im_bin_none - im_bin))/length(im_bin);

%% Recover image (no channel estimation)
% rec_im = reshape(rec_im_bin,9,numel(rec_im_bin)/8);
rec_im_none = reshape(rec_im_bin_none,8,numel(rec_im_bin_none)/8);
rec_im_none = uint8(bin2dec(rec_im_none'));
rec_im_none = reshape(rec_im_none, size(im));
end
