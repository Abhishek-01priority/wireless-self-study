%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% File to simulate clarkes model. The code has been written 
% as described in "Wireless Communications - Rappaort" page 183
% the transmitted signal should be convolved with 'r' to generate
% mulitpath scenario
% Input args: N -> Number of frequency points in doppler spectrum
%             fm -> maximum doppler shift
% Learnings:
%  1) the doppler spectrum reaches infinity at f = +- fm. In order
%     to avoid that, spectrum at f = +-fm is calculated by taking 
%     slope.
%  2) clarke model gives complex channel gains as output that the 
%     signal should be modulated with.
%  3) After coding I found that `r` has very low variance which implied
%     signal apart from deep fade events was going to be serverly faded.
%     I also found a video [Lecture 29 Introduction to Wireless and
%     Cellular system], where the variance of `r` is supposed to be 
%     0.5. Hence, `r` is multiplied by a gain that makes its variance
%     0.5.
%  4) In clarke's model, the fading channel is created for a duration 
%     of Nifft * Ts seconds. This duration is dependent on max doppler
%     as well as number frequency points in doppler spectrum. Hence, 
%     to create a fading for only 1ms, a huge number of Nifft is required.
%     And to generate huge amount of Nifft, tweaking wrt N and fm is needed.
%     This interdependcies are a major drawback of clarke's model.
% =========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  
  r_squared = ifft_I_noise + ifft_Q_noise;
  
  r = sqrt(r_squared);
  
  r = 1 / (sqrt ( 2 * var(real(r)) )) * real(r) + ...
      1 / (sqrt ( 2 * var(imag(r)) )) * imag(r);

end% function end
