% EE-442 Wireless Receivers: algorithms and architectures
% Final Project: OFDM Audio Transmission System
% EPFL Fall 2022
% Authors: Dong Chu, Han Miao, Huan-Ying Yeh

%% Visualization_05_Training_Intervals
% Change the interval of pilot/training symbol insertion and compare the
% resulting experimental BER with the audio transmission system

% Plotted experimental data.
x = [1 2 3];

ber_block = [0.0432, 0.2931, 0.4707];
ber_comb = [0.0013, 0.0021, 0.0076];

% Plot on the axis
figure();
hold on 
plot(x, ber_block , 'o-',  'Color','#0d9e20','LineWidth', 2.5);
plot(x, ber_comb ,  'o-',  'Color', '#c400ac','LineWidth', 2.5);

xlim([0 4]);
ylim([0 1]);

set(gca, 'YScale', 'log')
ticks = sort([ber_block, ber_comb]);
yticks(ticks)

xticks(x);
xticklabels({'1:1','1:2','1:4'})

xlabel('Train-Data Ratio');
ylabel('BER')
%ytickformat('percentage')

legend('Block', 'Comb', 'Location', 'Best');
grid on

title('Microphone Movement vs. BER Experimental Values');
saveas(gcf, 'plots/task03/Train_Intervals_vs_BER.png');
