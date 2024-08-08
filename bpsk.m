%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% This file does basic bpsk modulation and demodulation
% The implementation however was NOT SO BASIC. I struggled a lot to get 
% this to work.
% =========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

%#### CONFIGURATIONS ####%
% (fs) sampling frequency can be any aribitiary value
fs              = 1e3;
% (fc) a blog mentioned carrier frequency should be 1/4 of sampling frequency.
% fc is higher than symbol frequency (fs/16).
fc              = fs/4;
% (nspan) is the number of symbols the SRRC filter should span.
% nspan = 4 implies that symbols X, X+1, X+2, X+3 will be under the SRRC. X+1 and X+2 will
% come under SRRC major lobe where as X and X+3 in side lobes
nspan           = 4;
% (num_syms) is the total number of symbols for which the data should be generated.
% Ideally, the total number of bits for num_syms should be mod-order x num_syms and since
% bpsk has mod-order of 1, it doesn't really matter.
num_syms        = 4;
% (roll_off) is the SRRC parameter that governs the passband to stopband transition.
% the total bandwidth of system is (1 + roll_off)fs/samples_per_sym.
roll_off        = 0.25;
% (samples_per_sym) is the upsampling factor.
samples_per_sym = 16;

% carrier is generated for the upsampled sequence. This actually is the sampled version of carrier.
carrier = cos ( 2*pi*fc* ( 0 : samples_per_sym-1 ) / samples_per_sym );
fullcarrier = repmat(carrier, 1, num_syms);

%#### BPSK TRANSMITTER ####%
bits = randi([0 1], [1 num_syms]);
input_bits = [bits zeros(1, nspan/2)]; % to flush out the filter delay
bpsk_syms = 2*input_bits - 1;

% upsample the bpsk_syms to samples_per_sym symbols by adding 0s between consecutive symbols
zero_mat = zeros(samples_per_sym-1, num_syms + nspan/2);
bpsk_syms = [bpsk_syms; zero_mat];
bpsk_syms = reshape(bpsk_syms, 1, []);

% Square Root raised cosine pulse shaping
%   number of taps = nspan * samples_per_sym + 1. +1 is for middle peak tap.
%   Initially I wanted to derive SRRC from RRC by taking square root of RRC spectrum.
%   But, the resultant spectrum was incorrect. There are also no proper resources for it
%   so I dropped the idea and directly took the time domain formula for SRRC.
t = -nspan/2 : 1/samples_per_sym : nspan/2; % length(t) = total number of taps
%rrc = sinc(t) .* cos(roll_off * pi * t) ./ (1 - (2 * roll_off * t).^2 );
%rrc( find ( abs(rrc) == Inf ) ) = 0;
srrc = ( 2 * roll_off./(pi * sqrt(samples_per_sym)) ) .* ( cos( (1 + roll_off)*pi*t ) + ...
          sin ( (1 - roll_off)*pi*t ) ./ (4 * roll_off .* t) ) ./ ( 1 - (4 * roll_off * t).^2 );
% abs(srrc) will be Inf at zero crossings. Hence making those points 0
srrc( find (abs(srrc) == Inf) ) = 0;
% At middle peak tap, srrc = NaN because at that point t = 0. Therefore, the value at 
% middle peak tap is calculated by equation of line.
peak_loc = nspan/2 * samples_per_sym + 1;
m = ( srrc(peak_loc-1) - srrc(peak_loc-2) ) / ( t(peak_loc-1) - t(peak_loc-2) );
srrc(peak_loc) = m * ( t(peak_loc) - t(peak_loc-2) ) + srrc(peak_loc-2);

% convolving to make bpsk symbols continuous
%   The pulse shaping filter must be convolved with input signal.
%   This type of filter is called TX filter. But "convolving" doesn't necessarily
%   mean using conv() function. conv() perfoms linear convolution and output length is
%   length(filter) + length(input signal) - 1. Where as filter() output length is
%   length(input signal). For this reason, using filter was more intutive for me. 
%   However, filter() adds a filter delay that needs to be taken into consideration. 
%   If you see stem(srrc), the peak is at npsan/2*samples_per_sym+1 location and not at 1.
%   Implying a filter delay of nspan/2*samples_per_sym. It is very very important to 
%   remember that this filter delay should always be compensated whenever using filter()
pulse_shaped_bpsk = filter(srrc, 1, bpsk_syms);
fltDelay = (nspan / 2) * samples_per_sym;
pulse_shaped_bpsk = pulse_shaped_bpsk(fltDelay + 1: end);

% modulating with carrier frequency
mod_pulse_bpsk = pulse_shaped_bpsk .* fullcarrier;

%#### BPSK RECEIVER ####%

% De-Modulating BPSK
demod_pulse_bpsk = mod_pulse_bpsk .* fullcarrier;

% convolving received signal with Square root raised cosine
%   After demodulation, demod_pulse_bpsk will have a frequency component 2fc.
%   This has to be removed before next processing implying it has to be passed through
%   low pass filter. However, SRRC itself has a low pass nature and hence, additional
%   lpf is not needed here.
%   We perform match filtering by convolving received signal with SRRC and to compensate
%   for filter delay demod_pulse_bpsk is appended with zeros. If you notice, it is appended
%   with nspan/2*samples_per_sym and not npsan. This is because demod_pulse_bpsk is an upsampled 
%   sequence that is upsampled by samples_per_sym.
demod_pulse_bpsk = [demod_pulse_bpsk zeros(1, nspan/2 * samples_per_sym)];
matchFilteredBpsk = filter(srrc, 1, demod_pulse_bpsk);
matchFilteredBpsk = matchFilteredBpsk(nspan/2 * samples_per_sym + 1: end);

% sampling
outbits = matchFilteredBpsk(1 : samples_per_sym : end);
outbits(find(outbits > 0)) = 1;
outbits(find(outbits <= 0)) = 0;

bits
outbits

figure;
plot(pulse_shaped_bpsk,'-ob' ,'MarkerFaceColor', 'b')
hold on;
plot(matchFilteredBpsk,'-^r', 'MarkerFaceColor', 'r')
legend('pulse shaped bpsk', 'Match filtered bpsk')
grid minor

figure;
N = 1024;
plot(fs/N * (-N/2:N/2-1), dbaf(pulse_shaped_bpsk, N))
hold on
plot(fs/N * (-N/2:N/2-1), dbaf(matchFilteredBpsk, N), 'r')
legend('pulse shaped bpsk', 'Match filtered bpsk')
grid minor
