% EE-442 Wireless Receivers: algorithms and architectures
% Final Project: OFDM Audio Transmission System
% December 2022
% Authors: Dong Chu, Han Miao, Huan-Ying Yeh

% Audio Experiments: Comb Type Digital Transmitter
function [txsignal conf] = tx(txbits,conf,k)
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete transmitter
%   consisting of:
%       - modulator
%       - pulse shaping filter
%       - up converter
%   in digital domain.
%
%   Parameters:
%   txbits  : Information bits
%   conf    : Universal configuration structure
%   k       : Frame index
%
%   Outputs:
%   txsignal
%   conf
%
%% Parameters
n_sc = conf.sc; % number of sub-carriers
os_factor = conf.os_factor; % oversampling factor
L = conf.L; % length of cyclic prefix
os_trsym = conf.os_trsym; % insert training symbols every (os_trsym-1) data symbols

%% Generate OFDM symbols (training + data)
% QPSK modulation for OFDM symbols
txbits = reshape(txbits,[],2);
GreyMap_QPSK = 1/sqrt(2) * [-1-1j, -1+1j, 1-1j, 1+1j];
conf.GreyMap_QPSK = GreyMap_QPSK;
txsymbols = GreyMap_QPSK(bi2de(txbits,'left-msb')+1);

% known training symbols
num_train_symbol = ceil(n_sc/os_trsym); % number of training symbol datapoint every OFDM symbols
train_symbols = randi([0,1], 1, num_train_symbol); % N/2 (128) training symbol
train_symbols = 1 - 2*train_symbols; % BPSK
conf.tr_syms = train_symbols; % saving BPSK training symbols

% data symbols
num_data_symbol = n_sc - num_train_symbol; % nymber of data symbol datapoint every OFDM symbols
conf.num_data_symbol = num_data_symbol;
conf.nsyms = ceil(length(txsymbols)/num_data_symbol); % number of complete OFDM symbols, inserted noise at the last symbol
conf.extra_pnts = conf.nsyms * num_data_symbol - length(txsymbols); % number of padded random QPSK symbols at the tail
data_symbols = [txsymbols, GreyMap_QPSK(randi([1,4],1, conf.extra_pnts))]; % pad the last incomplete OFDM symbol

% paralize training and data symbols
train_symbols = repmat(train_symbols.', 1, conf.nsyms); % repetition for of training symbols
data_symbols = reshape(data_symbols, num_data_symbol, conf.nsyms);

% concatenate training and data OFDM symbols
freq_OFDMsymbols = zeros(n_sc, conf.nsyms);
idx_all = 1:n_sc; 
idx_train = 1:os_trsym:n_sc; % training symbols from the 1st position
idx_data = idx_all(~ismember(idx_all, idx_train));
conf.idx_train = idx_train;
conf.idx_data = idx_data;

for i = 1:conf.nsyms
    freq_OFDMsymbols(idx_train,i) = train_symbols(:,i);
    freq_OFDMsymbols(idx_data,i) = data_symbols(:,i);
end

%% Generate time-domain symbols
% Init
time_OFDMsymbols = zeros((n_sc * os_factor + L) * conf.nsyms, 1);

% IFFT and append cyclic prefix
for i = 0:conf.nsyms-1
    conf.beforeifft = freq_OFDMsymbols(:,i+1);
    IFFTsymbol = osifft(freq_OFDMsymbols(:,i+1),os_factor).';
    prefix = IFFTsymbol(end-L+1:end);
    time_OFDMsymbols(1 + i*(n_sc*os_factor+L) : (i+1)*(n_sc*os_factor+L)) = [prefix, IFFTsymbol];  
end

%% Generate the preamble
conf.frame_sync_length = 100;
preamble = 1 - 2*lfsr_framesync(conf.frame_sync_length); % generate in BPSK
conf.preamble = preamble;
up_preamble = upsample(preamble,os_factor); % upsampled preamble

%% Pulse shaping
% This is to ensure that the pulse shape will have zero value at next block
rrc_len = os_factor - 1;
rolloff = 0.22;
pulse = rrc(os_factor,rolloff, rrc_len);

% Convolving the upsampled sequence with the pulse shape filter
preamble = conv(up_preamble, pulse.','same');

% Power normalization
pwr_OFDM = norm(time_OFDMsymbols,2)/sqrt(length(time_OFDMsymbols));
pwr_preamble = norm(preamble,2)/sqrt(length(preamble));
preamble = preamble / pwr_preamble;
time_OFDMsymbols = time_OFDMsymbols / pwr_OFDM;

%% Concatenate preamble and OFDM symbols
OFDM_symbols = [preamble; time_OFDMsymbols];

%% Up-conversion
time = 1/conf.f_s:1/conf.f_s:length(OFDM_symbols)/conf.f_s; % same length as OFDM symbols
txsignal = OFDM_symbols.*exp(1j*2*pi*conf.f_c*time.');

txsignal = real(txsignal);