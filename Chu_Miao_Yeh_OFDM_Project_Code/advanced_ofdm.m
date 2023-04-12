% EE-442 Wireless Receivers: algorithms and architectures
% Final Project: OFDM Audio Transmission System
% EPFL Fall 2022
% Authors: Dong Chu, Han Miao, Huan-Ying Yeh

%% advanced_ofdm.m
% Simulation Helper Function
% Transmits and receives data with the specified training type and reyleigh noise channel
function [txbits, rxbits, rxsymbols, rx_img, ber] = advanced_ofdm(im, training_type, snr, n_taps, ray_channel)
%
% Calls the tx() and rx() functions of our EE-442 OFDM receivers (QPSK modulation)
% Apply AWGN and flat fading of desired SNR and #taps
%
% Params:
% im - input image file
% training_type - 'block' or 'comb' or 'viterbi' or 'vit'
% snr - SNR
% n_taps - fading tap numbers
% ray_channel - [Optional] the ray_channel object
%
% Outputs:
% - txbits: transmitted bits from the image
% - rxbits: received bits
% - rxsymbols: received symbols
% - rx_img: the received image
% - ber: bit error rate

%% Load raw data: image into txbits
bits = reshape((dec2bin(typecast(im(:), 'uint8'), 8) - '0').', 1, []);
txbits = bits';

%% Set up the config
conf.f_s     = 48000;   % sampling rate
conf.os_trsym = 2; % [input] insert training symbols every (os_trsym-1) data symbols
conf.f_sym   = 100;     % symbol rate
conf.nframes = 1;       % number of frames to transmit
conf.nbits   = size(txbits, 1);    % [input] number of bits
conf.modulation_order = 2; % Use QPSK
conf.f_c     = 8000; % carrier frequency
conf.sc = 256; % [input] number of subcarriers. Default = 256
conf.spacing = 4.6875;

% Spacing
conf.os_factor = conf.f_s/conf.sc/conf.spacing; % os_factor = 40
conf.L = conf.sc*conf.os_factor/2; % the length of cyclic prefix is half of the one of OFDM symbols
conf.L_h = 6; % [input] channel length
conf.framelen = 2; % [input] the number of OFDM symbols in a frame counting the training symbol

conf.npreamble  = 100;  % [Input] this could change
conf.bitsps     = 16;   % bits per audio sample
conf.offset     = 0;

%% Call the tx() based on pilot training type [Add Viterbi later]
switch(training_type)
    case 'block'     
        [~, txsignal_complex, ~, conf] = block_tx(txbits,conf,conf.nframes);
    case 'comb'   
        [~, txsignal_complex, ~, conf] = comb_tx(txbits,conf,conf.nframes);
    case 'viterbi'   
        conf.viterbi = 1;
        [~, txsignal_complex, ~, conf] = viterbi_tx(txbits,conf,conf.nframes);
end
    
%% [This section replaces the actual audio transmission & recording process]
% Add AWGN
x_s = txsignal_complex;
data_pwr = mean(abs(x_s.^2));

% Add noise
noise_pwr = data_pwr/10^(snr/10);
noise = normrnd(0,sqrt(noise_pwr/2),size(x_s))+normrnd(0,sqrt(noise_pwr/2),size(x_s))*1i;
x_s_noise = x_s + noise;

% Apply fading channel
g = exp(-(0:n_taps-1));
g = g/norm(g);
x_s_noise_fading = conv(x_s_noise,g,'same');

% Deal with Rayleigh Channels if applicable
if isa(ray_channel, 'comm.RayleighChannel') % Then a ray channel is supplied
    x_s_noise_fading = ray_channel(x_s_noise);
end

% Set up the array for transmission
normtxsignal = real(x_s_noise_fading);
rawtxsignal = [ zeros(conf.f_s,1) ; normtxsignal ;  zeros(conf.f_s,1) ]; % add padding before and after the signal
rawtxsignal = [  rawtxsignal  zeros(size(rawtxsignal)) ]; % add second channel: no signal

% [Audio transmits noisily]
rxsignal = rawtxsignal(:,1);

%% The rx() part
% Create struct to store results
res.biterrors   = zeros(conf.nframes,1);
res.rxnbits     = zeros(conf.nframes,1);

% Call rx() to recover the time-domain bits

%% Call the rx() to recover t-d bits basef on the training type
switch(training_type)
    case 'block'     
        [rxsymbols, rxbits, ~] = block_rx(rxsignal,conf);
    case 'comb'   
        [rxsymbols, rxbits, ~] = comb_rx(rxsignal,conf);
    case 'viterbi'   
        [rxsymbols, rxbits, ~] = viterbi_rx(rxsignal,conf);
end

%% Convert the received bits into image
orig_class = class(im);
orig_size = size(im);
rx_img = reshape(typecast(uint8(bin2dec(char(reshape(rxbits', 8, [])+'0').')), orig_class), orig_size);

%% Store BER
k = 1;
res.rxnbits(k)      = length(rxbits);
res.biterrors(k)    = sum(rxbits ~= txbits);
ber = sum(res.biterrors)/sum(res.rxnbits);

end
