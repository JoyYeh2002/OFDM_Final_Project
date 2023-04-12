% EE-442 Wireless Receivers: algorithms and architectures
% Final Project: OFDM Audio Transmission System
% EPFL Fall 2022
% Authors: Dong Chu, Han Miao, Huan-Ying Yeh

%% Visualization_04_Experiment_Plots
% Plot experimental data collected by physically moving the speakers back
% and forth and transmitting the audio data

x = [1 2 3];

ber_block = [0.0003, 0.0432, 0.3215];
ber_comb = [0.0000001, 0.0013, 0.0549];

% Plot on the axis
figure();
hold on 
plot(x, ber_block, 'o-',  'Color','#0d9e20','LineWidth', 2.5);
plot(x, ber_comb,  'o-',  'Color', '#c400ac','LineWidth', 2.5);

xlim([0 4]);
ylim([0 1]);

set(gca, 'YScale', 'log')
ticks = sort([ber_block, ber_comb]);
yticks(ticks)

xticks(x);
xticklabels({'No Movement','Slight Movement','Strong Movement'})

ylabel('BER')
% ytickformat('percentage')

legend('Block', 'Comb', 'Viterbi', 'Location', 'Best');
grid on

title('Microphone Movement vs. BER Experimental Values');
saveas(gcf, 'plots/task03/Exp_Movement_vs_BER.png');
