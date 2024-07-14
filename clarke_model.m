%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==================================================================
% File to simulate clarkes model. The code has been written 
% as described in "Wireless Communications - Rappaort" page 183
% the transmitted signal should be convolved with 'r' to generate
% mulitpath scenario
% Input args: N -> Number of FFT/IFFT bins
%             fm -> maximum doppler shift
% ==================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r, Nifft] = clarke_model(N, fm)
  Ts = 5.15e-6; % sampling time
  N = N; % number of frequency domain points to represent baseband doppler
  fm = fm; % (in Hz) maximum doppler shift
  fc = 0; % carrier frequency
  f = linspace(fc-fm, fc+fm, N); % frequency points are spread between maximum doppler shift

  delta_f = 2 * fm / (N - 1);
  T = 1/delta_f; % time duration of a fading waveform
  
  % generating complex gaussian noise
  % conjugate noise components lie in negative frequency
  noise_I_positive_freq = randn(1,N/2) + 1i * randn(1,N/2); % complex gaussian
  noise_I = [conj(noise_I_positive_freq) noise_I_positive_freq];
  
  noise_Q_positive_freq = randn(1,N/2) + 1i * randn(1,N/2); % complex gaussian
  noise_Q = [conj(noise_Q_positive_freq) noise_Q_positive_freq];
  
  % fading spectrum
  S_baseband = zeros(size(f));
  S_baseband = 1.5 ./ (pi * fm * sqrt(1 - ((f - fc) ./ fm).^2));
  S_baseband(1) = (S_baseband(3) - S_baseband(2))/ (f(3) - f(2)) * (f(1) - f(2)) + S_baseband(2);
  S_baseband(end) = (S_baseband(end-2) - S_baseband(end-1))/ (f(end-2) - f(end-1)) * (f(end) - f(end-1)) + S_baseband(end-1);
  
  %%% at this point all pre-requisties are ready
  
  doppler_spread_I_noise = sqrt(S_baseband) .* noise_I;
  doppler_spread_Q_noise = sqrt(S_baseband) .* noise_Q;
  
  Nifft = 2.^ ceil( log2(1 / (Ts * delta_f)));

  ifft_I_noise = ( ifft(doppler_spread_I_noise, Nifft) ) .^ 2;
  ifft_Q_noise = ( 1i * ifft(doppler_spread_Q_noise, Nifft) ) .^ 2;
  
  r_squared = (sqrt( 1 / (2*var(ifft_I_noise)) ) * ifft_I_noise) + ...
              (sqrt( 1 / (2*var(ifft_Q_noise)) ) * ifft_Q_noise);
  
  r = sqrt(r_squared);
end% function end
