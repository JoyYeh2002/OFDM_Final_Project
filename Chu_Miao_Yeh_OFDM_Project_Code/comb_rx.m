% EE-442 Wireless Receivers: algorithms and architectures
% Final Project: OFDM Audio Transmission System
% December 2022
% Authors: Dong Chu, Han Miao, Huan-Ying Yeh

%% Rayleigh Channel Simluations: Comb Type Digital Receiver
function [freq_OFDM_symbols, rxbits, conf] = comb_rx(rxsignal,conf,k)
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

% % dummy 
% rxbits = zeros(conf.nbits,1);

nsyms = conf.nsyms; % number of received OFDM symbols, half of which are training symbols.
L = conf.L; % length of cyclic prefix
n_sc = conf.sc; % number of subcarrier
os_factor = conf.os_factor; % oversampling factor
L_h = conf.L_h; % channel length

%% down conversion
time = 1/conf.f_s : 1/conf.f_s : length(rxsignal)/conf.f_s; % same length as OFDM symbols
rxsignal = rxsignal.*exp(-1j*2*pi*conf.f_c*time.');

% [!] Plot these later
% figure()
% RX = fft(rxsignal);
% t = (0:length(rxsignal)-1) / length(rxsignal)* conf.f_s;
% plot(t,abs(RX))
% title('FFT Mgnitude of the RX Signal');

%% lowpass 
f = conf.spacing * n_sc * 3;% corner frequency

[rxsignal] = ofdmlowpass(rxsignal,conf,f);

%% frame synchronization
[str_idx] = frame_sync(rxsignal,conf); % starting idx = 48048
% str_idx = str_idx+2;
rx_symbols = rxsignal(str_idx : str_idx + nsyms * (conf.os_factor * n_sc+L) - 1);

%% remove cyclic prefix and convert back to freq domain via FFT
OFDM_symbols = zeros(n_sc, nsyms );
for i = 0 : nsyms-1
    OFDM_symbol = rx_symbols(i * (n_sc * os_factor + L) + L + 1: (i+1) * (n_sc * os_factor + L)); % one OFDM symbol
%     if i==0
%         conf.after = OFDM_symbol; % for debugging
%     end
    OFDM_symbol = osfft(OFDM_symbol, os_factor); % downsampled during osFFT
    OFDM_symbols(:,i+1) = OFDM_symbol;    
end

%% channel estimation and removing trainign symbols
reverted_OFDM_symbols = zeros(n_sc, nsyms);
freq_OFDM_symbols = zeros(conf.num_data_symbol*nsyms,1);

FFT_mat = dftmtx(n_sc);
F = FFT_mat(conf.idx_train,1:L_h); % N/2*L FFT matrix
T = diag(conf.tr_syms); 

for i = 1:nsyms
    y = OFDM_symbols(conf.idx_train, i);
    h_hat = ((F') * (T') * T * F)\ ((F') * T' * y);
    H_hat = fft(h_hat, n_sc);
% 
%     figure()
%     plot(1:length(H_hat),abs(H_hat))
%     title("The channel of the" + num2str(i) + "-th symbol");
%     
    reverted_OFDM_symbols(:,i) = OFDM_symbols(:,i) ./H_hat;
    freq_OFDM_symbols((i-1)*conf.num_data_symbol+1 : i*conf.num_data_symbol) = reverted_OFDM_symbols(conf.idx_data,i);
end

%% Demodulation
freq_OFDM_symbols = freq_OFDM_symbols(1: end-conf.extra_pnts); % remove inserted extra points
GreyMap_QPSK = conf.GreyMap_QPSK;
[~,idx] = min( repmat(GreyMap_QPSK,length(freq_OFDM_symbols),1) - repmat(freq_OFDM_symbols, 1, 4), [], 2 );
rxbits = de2bi(idx-1, 'left-msb');
rxbits = rxbits(:);

