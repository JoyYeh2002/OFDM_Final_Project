% EE-442 Wireless Receivers: algorithms and architectures
% Final Project: OFDM Audio Transmission System
% EPFL Fall 2022
% Authors: Dong Chu, Han Miao, Huan-Ying Yeh
%
%% Speaker and Microphone System Measurements
%% Main Driver Script
%
% Workflow:
% 1. Load, rescale, and transmit a 44x51-pixal gray-scale image of a flower
% 2. Use the 3 operating modes to transmit the bits data with audio
% 3. Pass input data into tx() and rx(), then calculate PER and BER
%
% Audio Recording
%   3 operating modes:
%   - 'matlab' : generic MATLAB audio routines (unreliable under Linux)
%   - 'native' : OS native audio system
%       - ALSA audio tools, most Linux distrubtions
%       - builtin WAV tools on Windows
%   - 'bypass' : no audio transmission, takes txsignal as received signal

clc;clear;close all;

% Configuration Values
conf.audiosystem = 'matlab'; % Values: 'matlab','native','bypass'

conf.f_s     = 48000;   % sampling rate  
conf.os_trsym = 5; % [input]  insert training symbols every (os_trsym-1) data symbols
conf.f_sym   = 100;     % symbol rate
conf.nframes = 1;       % number of frames to transmit
conf.nbits   = 8000;    % [input]  number of bits 
conf.modulation_order = 2; % BPSK:1, QPSK:2
conf.f_c     = 8000; % carrier frequency
conf.sc = 256; % number of subcarriers
conf.spacing = 4.6875; 

% spacing 
conf.os_factor = conf.f_s/conf.sc/conf.spacing; % os_factor = 40
conf.L = conf.sc*conf.os_factor/2; % the length of cyclic prefix is half of the one of OFDM symbols
conf.L_h = 6; % [input] channel length

conf.npreamble  = 100;
conf.bitsps     = 16; % bits per audio sample
conf.offset     = 0;
conf.showChannel = 0; % whether to visualize the estimated channel

% Init Section
% all calculations that you only have to do once
% conf.os_factor  = conf.f_s/conf.f_sym;
% if mod(conf.os_factor,1) ~= 0
%    disp('WARNING: Sampling rate must be a multiple of the symbol rate'); 
% end
% conf.nsyms      = ceil(conf.nbits/conf.modulation_order);

% Initialize result structure with zero
res.biterrors   = zeros(conf.nframes,1);
res.rxnbits     = zeros(conf.nframes,1);

% Results
for k=1:conf.nframes
     
    % Load input image
    img_name = 'data\flower.png';
    
    resize_scale = 0.1;
    I = imread(img_name);
    I = imresize(I, resize_scale);
    I = rgb2gray(I);
    
    figure()
    imshow(I); title('Input Image');
    
    % Transfom image into input bits
    bits = reshape((dec2bin(typecast(I(:), 'uint8'), 8) - '0').', 1, []);
    txbits = bits';

    % For transmission experiments with audio equipment, call the comb tx() transmit function
    [txsignal conf] = tx(txbits,conf,k);
    
    % % % % % % % % % % % %
    % Audio Transmission Begins
    
    % normalize values
    peakvalue       = max(abs(txsignal));
    normtxsignal    = txsignal / (peakvalue + 0.3);
    normtxsignal = txsignal;
    
    % create vector for transmission
    rawtxsignal = [ zeros(conf.f_s,1) ; normtxsignal ;  zeros(conf.f_s,1) ]; % add padding before and after the signal
    rawtxsignal = [  rawtxsignal  zeros(size(rawtxsignal)) ]; % add second channel: no signal
    txdur       = length(rawtxsignal)/conf.f_s; % calculate length of transmitted signal
   
    % Platform native audio mode 
    if strcmp(conf.audiosystem,'native')
        
        % Windows WAV mode 
        if ispc()
            disp('Windows WAV');
            wavplay(rawtxsignal,conf.f_s,'async');
            disp('Recording in Progress');
            rawrxsignal = wavrecord((txdur+1)*conf.f_s,conf.f_s);
            disp('Recording complete')
            rxsignal = rawrxsignal(1:end,1);

        % ALSA WAV mode 
        elseif isunix()
            disp('Linux ALSA');
            cmd = sprintf('arecord -c 2 -r %d -f s16_le  -d %d in.wav &',conf.f_s,ceil(txdur)+1);
            system(cmd); 
            disp('Recording in Progress');
            system('aplay  out.wav')
            pause(2);
            disp('Recording complete')
            rawrxsignal = audioread('in.wav');
            rxsignal    = rawrxsignal(1:end,1);
        end
        
    % MATLAB audio mode
    elseif strcmp(conf.audiosystem,'matlab')
        disp('MATLAB generic');
        playobj = audioplayer(rawtxsignal,conf.f_s,conf.bitsps);
        recobj  = audiorecorder(conf.f_s,conf.bitsps,1);
        record(recobj);
        disp('Recording in Progress');
        playblocking(playobj)
        pause(0.5);
        stop(recobj);
        disp('Recording complete')
        rawrxsignal  = getaudiodata(recobj,'int16');
        rxsignal     = double(rawrxsignal(1:end))/double(intmax('int16')) ;
        
    % Bypass mode doesn't pass the signal through anything. Returns the
    % same output as input
    elseif strcmp(conf.audiosystem,'bypass')
        rawrxsignal = rawtxsignal(:,1);
        rxsignal    = rawrxsignal;
    end

    % Audio Transmission Ends
    % % % % % % % % % % % %
    
    % For audio equipment experiments, call the comb rx() receiver function
    [rxbits, conf]       = rx(rxsignal,conf);

    % Convert received bits into image
    orig_class = class(I);
    orig_size = size(I);
    im_reconstructed = reshape(typecast(uint8(bin2dec(char(reshape(rxbits', 8, [])+'0').')), orig_class), orig_size);
    
    figure()
    imshow(im_reconstructed); title('Received Image');
    
    res.rxnbits(k)      = length(rxbits);  
    res.biterrors(k)    = sum(rxbits ~= txbits);
    
end

% Display the PER and BER
per = sum(res.biterrors > 0)/conf.nframes
ber = sum(res.biterrors)/sum(res.rxnbits)