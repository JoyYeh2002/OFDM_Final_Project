% EE-442 Wireless Receivers: algorithms and architectures
% Final Project: OFDM Audio Transmission System
% December 2022
% Authors: Dong Chu, Han Miao, Huan-Ying Yeh

%% [Appendix] Viterbi Digital Receiver
function [data_symbols, rxbits, conf] = viterbi_rx(rxsignal,conf,k)

% Digital Receiver
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete causal
%   receiver in digital domain.
%
%   rxsignal    : received signal
%   conf        : configuration structure
%   k           : frame index
%
%   Outputs
%
%   rxbits      : received bits
%   conf        : configuration structure
%
nsyms = conf.nsyms + 1; % a training symbol in the current version
L = conf.L; % length of cyclic prefix
n_sc = conf.sc; % number of subcarrier
os_factor = conf.os_factor;

%% down conversion
time = 1/conf.f_s : 1/conf.f_s : length(rxsignal)/conf.f_s; % same length as OFDM symbols
rxsignal = rxsignal.*exp(-1j*2*pi*conf.f_c*time.');
rxsignal = ofdmlowpass(rxsignal,conf,conf.f_c)*2;

%% frame synchronization
[str_idx, phase_offset] = viterbi_frame_sync(rxsignal,conf); % starting idx
rx_symbols = rxsignal(str_idx : str_idx + conf.numtotal * (conf.os_factor * n_sc+L) - 1);

%% remove cyclic prefix and convert back to freq domain via FFT
OFDM_symbols = zeros(n_sc, conf.numtotal);
for i = 0 : conf.numtotal-1
    OFDM_symbol = rx_symbols(i * (n_sc * os_factor+L) + L +1: (i+1) * (n_sc * os_factor+L)); % one OFDM symbol
    OFDM_symbol = osfft(OFDM_symbol, os_factor); % downsampled during osFFT
    OFDM_symbols(:,i+1) = OFDM_symbol;    
end

Ch_est = OFDM_symbols(:,conf.train_index)./repmat(conf.tr_syms,1,conf.numtrain);
data_symbols = zeros(conf.sc, conf.nsyms);

if(conf.viterbi)
    for i = 0:size(Ch_est,2)-1

        theta_hat = zeros(conf.sc, conf.framelen);
        theta_hat(:,1) = angle(Ch_est(:,i+1));

        for k = 1:(conf.framelen-1)
            
            deltaTheta = 1/4*angle(-OFDM_symbols(:,i*conf.framelen+1+k).^4) + pi/2*(-1:4);

            [theta, ~] = min(abs(deltaTheta - theta_hat(:,k)),[],2);
            
            theta_hat(:,k+1) = mod(0.01*theta + 0.99*theta_hat(:,k), 2*pi);
            
            data_symbols(:,i*(conf.framelen-1)+k) = OFDM_symbols(:,i*conf.framelen+1+k) ./ abs(Ch_est(:,i+1)).*exp(-1j * theta_hat(:,k+1));

            if(i == size(Ch_est,2)-1 && k == mod(conf.nsyms,conf.framelen-1))
                break;
            end

        end
    end

else 

    %% Channel Estimation
    Ch_est = repelem(Ch_est,1,conf.framelen-1);
    
    CE_short = Ch_est(:,1:size(OFDM_symbols(:,not(conf.train_index)),2));
    
    
    %% Channel Correction
    data_symbols = OFDM_symbols(:,not(conf.train_index))./CE_short;

end

%% Demodulation
truncated = data_symbols(1:end-conf.extra_pnts); % flatten
truncated = transpose(truncated);
GreyMap_QPSK = conf.GreyMap_QPSK;

[~,idx] = min( repmat(GreyMap_QPSK,length(truncated),1) - repmat(truncated, 1, 4), [], 2 );

rxbits = de2bi(idx-1, 'left-msb');
rxbits = rxbits(:);