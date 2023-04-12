function [txsymbols, txsignal_complex, txsignal, conf] = block_tx(txbits,conf,k)
% Digital Transmitter
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete transmitter
%   consisting of:
%       - modulator
%       - pulse shaping filter
%       - up converter
%   in digital domain.
%
%   txbits  : Information bits
%   conf    : Universal configuration structure
%   k       : Frame index
%

% if given txsignal, have to quantize to 8 bits per sample first

% % dummy 400Hz sinus generation
% time = 1:1/conf.f_s:4; % 4 seconds?
% txsignal = 0.3*sin(2*pi*400 * time.');

%% parameters
n_sc = conf.sc; % number of sub-carriers

os_factor = conf.os_factor; % oversampling factor

L = conf.L; % length of cyclic prefix


%% Generate OFDM symbols (training + data)
% QPSK modulation for OFDM symbols
txbits = reshape(txbits,[],2);


%%% converting to QPSK %%%
GreyMap_QPSK = 1/sqrt(2) * [-1-1j, -1+1j, 1-1j, 1+1j];

conf.GreyMap_QPSK = GreyMap_QPSK;

txsymbols = GreyMap_QPSK(bi2de(txbits,'left-msb')+1);
%%% converting to QPSK %%%



conf.nsyms = ceil(length(txsymbols)/n_sc); % number of data OFDM symbols needed to contain the signal

conf.extra_pnts = conf.nsyms * n_sc - length(txsymbols); % number of padded random QPSK symbols at the tail

txsymbols = [txsymbols, GreyMap_QPSK(randi([1,4],1,conf.extra_pnts))]; % pad the last incomplete OFDM symbol

OFDM_data_symbols = reshape(txsymbols,n_sc,conf.nsyms); % Shaping this into a matrix

% known training OFDM symbols
OFDM_train_symbols = randi([0,1], n_sc,1);
conf.tr_syms = 1 - 2*OFDM_train_symbols; % BPSK

% Total number of OFDM symbols needed including the training symbols. 
conf.numtotal = ceil(conf.nsyms/(conf.framelen-1))+conf.nsyms; 

conf.numtrain = conf.numtotal-conf.nsyms;

% concatenate training and data OFDM symbols
freq_OFDMsymbols = zeros(n_sc,conf.numtotal);

column_index = [1:1:conf.numtotal]; 

conf.train_index = mod(column_index,conf.framelen)==1;

freq_OFDMsymbols(:,conf.train_index) = repmat(conf.tr_syms,1,sum(conf.train_index==1)); 
freq_OFDMsymbols(:,not(conf.train_index)) = OFDM_data_symbols;

%% Generate time-domain symbols

% Init
time_OFDMsymbols = zeros((n_sc * os_factor + L) * size(freq_OFDMsymbols,2), 1);

% IFFT and append cyclic prefix
for i = 0:size(freq_OFDMsymbols,2)-1
    IFFTsymbol = osifft(freq_OFDMsymbols(:,i+1),os_factor).';
    prefix = IFFTsymbol(end-L+1:end);
    time_OFDMsymbols(1 + i*(n_sc*os_factor+L) : (i+1)*(n_sc*os_factor+L)) = [prefix, IFFTsymbol]; 

end

%% Generate the preamble
conf.frame_sync_length = 100;
preamble = 1 - 2*lfsr_framesync(conf.frame_sync_length); % generate in BPSK
conf.preamble = preamble;
up_preamble = upsample(preamble,os_factor); % upsampled preamble

% Pulse shaping
% This is to ensure that the pulse shape will have zero value at the next ?
rrc_len = os_factor - 1; % remark 1: remain zeros between 2 rrc?
% rrc_len = os_factor*3-1; % remark 1: remain zeros between 2 rrc?
rolloff = 0.22;
pulse = rrc(os_factor,rolloff, rrc_len);

% Convolving the upsampled sequence with the pulse shape filter
preamble = conv(up_preamble, pulse.','same');

% power normalization
pwr_OFDM = norm(time_OFDMsymbols,2)/sqrt(length(time_OFDMsymbols));

pwr_preamble = norm(preamble,2)/sqrt(length(preamble));

preamble = preamble / pwr_preamble;

time_OFDMsymbols = time_OFDMsymbols / pwr_OFDM;

%% Concatenate preamble and OFDM symbols
OFDM_symbols = [preamble; time_OFDMsymbols];

%% up-conversion
time = 1/conf.f_s:1/conf.f_s:length(OFDM_symbols)/conf.f_s; % same length as OFDM symbols
txsignal_complex = OFDM_symbols.*exp(1j*2*pi*conf.f_c*time.');
txsignal = real(txsignal_complex);
end
