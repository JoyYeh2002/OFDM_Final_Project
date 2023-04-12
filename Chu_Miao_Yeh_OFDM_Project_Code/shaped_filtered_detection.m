function [filtered_rx_signal, pulse] = shaped_filtered_detection(rx_signal, rolloff, rx_filter_len, os_factor)
% input
% rx_signal: a column vector containing the received signal with awgn noise
% rolloff: roll-off factor of rrc filter, here we use 0.22
% rx_filter_len: length of the pulse for match filtering, here we use 12
% os_factor: the upsampling factor of system, here we use 4

% output: 
% filtered_rx_signal: a column vector
% pulse:  a column vector
%     filtered_rx_signal  = zeros(length(rx_signal)-2*rx_filter_len, 1);
    % define matched filter
     pulse = rrc(os_factor, rolloff, rx_filter_len);

    % filtering (one for loop)

     filtered_rx_signal = conv(rx_signal, pulse,'valid');


end