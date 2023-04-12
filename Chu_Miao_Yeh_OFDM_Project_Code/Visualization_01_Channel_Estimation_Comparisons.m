% EE-442 Wireless Receivers: algorithms and architectures
% Final Project: OFDM Audio Transmission System
% December 2022
% Authors: Dong Chu, Han Miao, Huan-Ying Yeh

%% Visualization Task 01: Channel Estimation Comparisons
% 
% Tunable inputs: 
% - input image
% - snr value and number of fading taps
% - various parameters of a Rayleigh fading channel
% 
% Outputs: 
% Plot a 5x2 subplot of the following:
% Input constellation, received constellation, with block chEst, comb
% chEst, and viterbi chEst
% The 2nd row has the corresponding image results 

clear all; close all; clc;

%% Initial set up
im = imread('data/mountain.png');
resize_scale = 0.2;
im = imresize(im, resize_scale);
im = rgb2gray(im);

%% Inputs: SNR, number of exponential fading taps, and Rayleigh channel
snr = 1; % [Input]
n_taps = 8; % [Input]

% Rayleigh Channel initialization
add_ray_channel = 'true'; % [Input} 'true', 'false';
sampleRate  = 48000;     % Sample rate of 20K Hz
maxDopplerShift  = 1; % Max Doppler shift of diffuse components (Hz)
delayVector = [1 2]*1e-8; % Discrete delays of four-path channel (s)
gainVector  = [-1 -1];  % Average path gains (dB)

% Create Rayleigh channel object if applicable
if strcmp(add_ray_channel, 'true')
    ray_channel = comm.RayleighChannel( ...
        'SampleRate',sampleRate, ...
        'PathDelays',delayVector, ...
        'AveragePathGains',gainVector, ...
        'NormalizePathGains',true, ...
        'MaximumDopplerShift',maxDopplerShift, ...
        'RandomStream','mt19937ar with seed', ...
        'Seed',10, ...
        'PathGainsOutputPort',true);
else
    ray_channel = 'false';
end

% Call the simple function
[X, X_hat_none, X_hat, rec_im_none, rec_im, ber_none, ber] = simple_ofdm(im, snr, n_taps, ray_channel);

% Calculate advanced_ofdm results
[block_txbits, block_rxbits, block_rxsymbols, block_rx_img, block_ber] = advanced_ofdm(im, 'block', snr, n_taps, ray_channel);
[comb_txbits, comb_rxbits, comb_rxsymbols, comb_rx_img, comb_ber] = advanced_ofdm(im, 'comb', snr, n_taps, ray_channel);
[vit_txbits, vit_rxbits, vit_rxsymbols, vit_rx_img, vit_ber] = advanced_ofdm(im, 'viterbi', snr, n_taps, ray_channel);

%% Generate plots
plot_figs = 1;
if plot_figs == 1
    
% Colors: 
plot_color = '#3892eb';
% 1: #7F167F dark purple
% 2: #CB1C8D magenta
% 3: '#F06292'
f = figure;
f.Position = [100 100 1200 450];
    
%% 1.1 transmit constellation
subplot(2,5,1);
plot(X,'x','linewidth',2, 'color', plot_color, 'markersize',10);
xlim([-2 2]);
ylim([-2 2]);
xlabel('In-Phase')
ylabel('Qudrature')

title(sprintf('\\bfTransmit Constellation\n\\rm SNR: %.2f dB. QPSK\n Channel Taps: %d', snr, n_taps));
grid on

%% 1.2 Raw Recovered constellation
subplot(2,5,2);
plot(X_hat_none(1:100:end),'x', 'color', plot_color, 'markersize',3);
xlim([-2 2]);
ylim([-2 2]);
xlabel('In-Phase')
ylabel('Qudrature')

title(sprintf('\\bf Raw Received\n\\rm'));
grid on

%% 1.3 Rec Constellation With Block Type
subplot(2,5,3);
plot(block_rxsymbols(1:100:end),'x', 'color', plot_color, 'markersize',3);
xlim([-2 2]);
ylim([-2 2]);
xlabel('In-Phase')
ylabel('Qudrature')

title(sprintf('\\bf Block Type Estimation\n\\rm'));
grid on

%% 1.4 Rec Constellation With Comb Type
subplot(2,5,4);
plot(comb_rxsymbols(1:100:end),'x', 'color', plot_color, 'markersize',3);
xlim([-2 2]);
ylim([-2 2]);
xlabel('In-Phase')
ylabel('Qudrature')

title(sprintf('\\bf Comb Type Estimation\n\\rm'));
grid on

%% 1.5 Rec Constellation With Viterbi Type
subplot(2,5,5);
plot(vit_rxsymbols(1:500:end)./100000,'x', 'color', plot_color, 'markersize',3);
%xlim([-2 2]);
%ylim([-2 2]);
xlabel('In-Phase')
ylabel('Qudrature')

title(sprintf('\\bf Viterbi Type Estimation\n\\rm'));
grid on

%% 2.1 Transmit image
subplot(2,5,6);
imshow(im);
title(sprintf('\\bfTransmit Image\n \\rm%d x %d Pixels', 103, 103));
grid on

%% 2.2 Recovered image without estimation
subplot(2,5,7);
imshow(rec_im_none);
title(sprintf('\\bf Raw\n\\rmBER: %.2g',ber_none));

%% 2.3 Recovered image Block Type
subplot(2,5,8);
imshow(block_rx_img);
title(sprintf('\\bf Block Type\n\\rmBER: %.2g', block_ber));

%% 2.4 Recovered image Comb Type
subplot(2,5,9);
imshow(comb_rx_img);
title(sprintf('\\bf Comb Type\n \\rmBER: %.2g', comb_ber));

%% 2.5 Recovered image Viterbi Type
subplot(2,5,10);
imshow(vit_rx_img);
title(sprintf('\\bf Viterbi Type\n \\rmBER: %.3f', vit_ber));
end

%% Save the plots
saveas(gcf, [pwd, '\plots\task04\Y2_Rayleigh_', add_ray_channel, '_SNR_', num2str(snr), '_taps_', num2str(n_taps), '.png']);