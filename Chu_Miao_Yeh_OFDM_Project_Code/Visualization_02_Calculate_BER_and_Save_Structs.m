% EE-442 Wireless Receivers: algorithms and architectures
% Final Project: OFDM Audio Transmission System
% EPFL Fall 2022
% Authors: Dong Chu, Han Miao, Huan-Ying Yeh

%% Visualization Task 02: Calculate_BER_and_Save_Structs.m
%% Rayleigh Parameters vs. BER
% Change various parameters of the Rayleigh fading and evaluate the BER
% Save the results to a struct.

clear all; close all; clc;

% These parameters are fixed
snr = 1;
n_taps = 8;

% Initial set up
im = imread('data/mountain.png');
resize_scale = 0.2;
im = imresize(im, resize_scale);
im = rgb2gray(im);

% Create empty struct
BER = {};

%% Generate various testing points for the Rayleigh channels and save the BER results to a struct
fieldName = 'delay'; % 'SNR', 'doppler', 'n_taps'
var = {[0 0], [1 2]*1e-5, [1 2]*1e-3, [1 2]*1e-2, [0.1, 0.2]};
for idx = 1 : size(var, 2)
    sampleRate = 48000;     % Sample rate of 20K Hz
    maxDopplerShift = 0; % Max Doppler shift of diffuse components (Hz)
    delayVector = var{idx}; % Discrete delays of four-path channel (s)
    gainVector  = [-16 -9];  % Average path gains (dB)
    
    % Create the channel objects
    ray_channel = comm.RayleighChannel( ...
        'SampleRate',sampleRate, ...
        'PathDelays',delayVector, ...
        'AveragePathGains',gainVector, ...
        'NormalizePathGains',true, ...
        'MaximumDopplerShift',maxDopplerShift, ...
        'RandomStream','mt19937ar with seed', ...
        'Seed',10, ...
        'PathGainsOutputPort',true);
    
    % Call the simple function
    %[X, X_hat_none, ~, rec_im_none, ~, ber_none, ~] = simple_ofdm(im, snr, n_taps, ray_channel);
    
    % Calculate advanced_ofdm results
    [~, ~, ~, ~, ber_block] = advanced_ofdm(im, 'block', snr, n_taps, ray_channel);
    [~, ~, ~, ~, ber_comb] = advanced_ofdm(im, 'comb', snr, n_taps, ray_channel);
    [~, ~, ~, ~, ber_vit] = advanced_ofdm(im, 'viterbi', snr, n_taps, ray_channel);
    
    % Store the BER values into the struct
    BER.(fieldName){idx} = var(idx);
    BER.block{idx} = ber_block;
    BER.comb{idx} = ber_comb;
    BER.vit{idx} = ber_vit;
end

BER.(fieldName) = [0, 1e-7, 1e-6, 1e-5, 1e-4];
BER.block = cell2mat(BER.block);
BER.comb = cell2mat(BER.comb);
BER.vit = cell2mat(BER.vit);

%% This script is run several times with different parameters/fields of the
% Rayleigh channel object. 
% The Visualization_02_BER_Plots.m script visualizes all the BER data.
save(['data/Ray_BER_', fieldName, '.mat'], 'BER');
