function [after] = lowpass(before,conf)
% LOWPASS lowpass filter
% Low pass filter for extracting the baseband signal 
%
%   before  : Unfiltered signal
%   conf    : Global configuration variable
%
%   after   : Filtered signal
%
% Note: This filter is very simple but should be decent for most 
% application. For very high symbol rates and/or low carrier frequencies
% it might need tweaking.
%

h_lp=design(fdesign.lowpass('N,Fc',20,conf.f_symb*2,conf.f_s),'all','MinPhase',true);
after = filter(h_lp,before);
