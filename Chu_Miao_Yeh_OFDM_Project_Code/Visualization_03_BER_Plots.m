% EE-442 Wireless Receivers: algorithms and architectures
% Final Project: OFDM Audio Transmission System
% EPFL Fall 2022
% Authors: Dong Chu, Han Miao, Huan-Ying Yeh

%% Visualization Task 03: BER_Plots.m
% Load structs saved from the "Visualization02.m" script and plot

close all;
%% Load Struct
fileName = 'delay04'; % 'delay', 'doppler', 'n_taps', 'snr', 'SNR_Flat'

load(['data/Ray_BER_', fileName, '.mat']);
fieldName = 'delay';
% Plot on the axis
figure();
hold on 

plot(BER.(fieldName), BER.block, 'o-', 'LineWidth', 2.5);
plot(BER.(fieldName), BER.comb, 'o-', 'LineWidth', 2.5);
%plot(BER.(fieldName), BER.vit, 'o-', 'LineWidth', 2.5);

%set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

xlim([0 1.2]);
%ylim([0 0.5]);

% ticks = sort([BER.block, BER.comb]);
%xticks(BER.(fieldName))

%% Change the x-labels

% xlabel('SNR in flat-fading channel');
% xlabel('Max Doppler Shift(Hz)')
% xlabel('Number of Exponential Fading Taps');
xlabel('SNR with Rayleigh Channel');
ylabel('BER')

%ytickformat('percentage')

legend('Block', 'Comb', 'Viterbi', 'Location', 'Best');
grid on

%% Change the titles
title('Rayleigh Doppler Shift vs. BER Plot: Zoomed-in');
% title('SNR with Rayleigh vs. BER Plot');
% title('SNR with Rayleigh vs. BER Plot');

%% Save the plots
saveas(gcf, ['plots/task03/Ray_', fieldName, '_vs_BER.png']);
    
