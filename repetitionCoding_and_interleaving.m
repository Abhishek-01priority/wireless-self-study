%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================================================================
% This file investigates the effect of interleaving and repetition 
% coding on a fading channel. The fading channel is generated from 
% clarke's model. I could have looped over different 
% repetitionCodingRate values but that would make the file messy.
% Learnings:
%  1) It is important to fix the randn seed for AWGN noise since 
%     Rayleigh fading is already constant.
%  2) As the repetitionCodingRate increases the number of new symbols
%     decreases and hence the amount of data transmitted in given 
%     time reduces. Essentially data rate decreases with increase in
%     repetitionCodingRate.
%  3) Repetition and/or Interleaving benefits can only be realized 
%     via coherent accumulation. With coherent accumulation, the 
%     coding gain also can be observed
%  4) From plot, the Pe at high SNRs is decided by fading events.
%     Interleaving has clear advantage in those instances.
% ======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixing the randn seed.
V = randn("seed");

% generate clarkes model
% 22 samples in doppler spectrum and 100Hz of doppler shift
[h, total_number_of_samples] = clarke_model(22, 100); 

N = total_number_of_samples;
% transmitted BPSK signal
signalTX = rand(1, N) - 0.5;
signalTX(find(signalTX <= 0)) = -1;
signalTX(find(signalTX > 0))  = +1;
sigPower = sum((signalTX.^2)/length(signalTX));

% snr
snr_in_dB = linspace(-15,15,10);

% this value should always be a power of 2
reptitionCodingRate = 4;

% the data rate will decrease by reptitionCodingRate
noSamplesForRepInterCoding = total_number_of_samples / reptitionCodingRate;
repIntSignaltx = rand(1,noSamplesForRepInterCoding) - 0.5;
repIntSignaltx(find(repIntSignaltx <= 0)) = -1;
repIntSignaltx(find(repIntSignaltx > 0))  = 1;

% repetition coding: 
% Each sample is repeated reptitionCodingRate times.
% So total size of transmitted signal is reptitionCodingRate*N
rep_signalTX = reshape( repmat(repIntSignaltx, reptitionCodingRate, 1), 1, [] );

% interleaving after repetition coding
% the repeated symbols are spread accross coherence time.
intlvd_signalTX = repmat(repIntSignaltx, 1, reptitionCodingRate);

% channel simulation for repetition and interleaving
for s = 1:length(snr_in_dB)
  
  randn("seed", V);
  % AWGN noise addition. transmited signal power for repetition and interleaving is
  % equal to transmited power of non-repeated-non-interleaved signal.
  noisePwr = 10 ^ (-snr_in_dB(s)/10) * sigPower;
  noise = sqrt(noisePwr/2) .* ((randn(1,N)) + 1i * (randn(1,N)));
  noise = noise - mean(noise);

  % Idea is to ensure that repetition coding and interleaving undergo same 
  % AWGN noise and fading
  rx_repeated_signal = h .* rep_signalTX  + noise;
  rx_interleaved_signal = h .* intlvd_signalTX + noise;

  % signal estimation
  y_tilda_repetition = conj(h) .* rx_repeated_signal ./ sqrt(conj(h) .* h);
  y_tilda_interleaved = conj(h) .* rx_interleaved_signal ./ sqrt(conj(h) .* h);

  % coherent accumulation
  coherent_y_rep = sum ( reshape( y_tilda_repetition, reptitionCodingRate, [] ) ); % row-wise sum
  coherent_y_int = sum ( reshape( y_tilda_interleaved, noSamplesForRepInterCoding, [] ) , 2)'; % columnwise sum

  % signal detection
  signalRX_repetition = real(coherent_y_rep);
  signalRX_repetition(find(signalRX_repetition <= 0)) = -1;
  signalRX_repetition(find(signalRX_repetition > 0)) = +1;

  signalRX_interleaved = real(coherent_y_int);
  signalRX_interleaved(find(signalRX_interleaved <= 0)) = -1;
  signalRX_interleaved(find(signalRX_interleaved > 0)) = +1;

  % Probability of Error calculation
  Pe_repetition(s) = length(find(signalRX_repetition ~= repIntSignaltx)) / N;
  Pe_interleaved(s) = length(find(signalRX_interleaved ~= repIntSignaltx)) / N;

end % for end

% channel sim for non-repetition and non-interleaving
for s = 1:length(snr_in_dB)
  
  randn("seed", V);
  % AWGN noise addition
  noisePwr = 10 ^ (-snr_in_dB(s)/10) * sigPower;
  noise = sqrt(noisePwr/2) .* ((randn(1,N)) + 1i * (randn(1,N)));

  %h = sqrt(1/2) * ( randn(1,N) + 1i * randn(1,N));
  %h = h - mean(h);
  received_signal = h .* signalTX + noise;

  % signal detection
  y_tilda = conj(h) .* received_signal ./ sqrt(conj(h) .* h);
  signalRX = real(y_tilda);
  signalRX(find(signalRX <= 0)) = -1;
  signalRX(find(signalRX > 0)) = +1;
  %signalRX = floor(signalRX + 0.5);

  % Probability of Error calculation
  Pe_Fading(s) = length(find(signalRX ~= signalTX)) / N;

end % for end

% plot
semilogy(snr_in_dB, Pe_repetition, '-ko', 'MarkerFaceColor', 'k')
hold on;
semilogy(snr_in_dB, Pe_interleaved, '-r^', 'MarkerFaceColor', 'r')
semilogy(snr_in_dB, Pe_Fading, '-bs', 'MarkerFaceColor', 'b')
title("Analysis of performance of repetition with/without interleaving in fading channel")
legend("Repetition coding", "Repetition and Interleaving", "No repetition and no interleaving")
xlabel("SNR in dB")
ylabel("Probability of Error")
grid minor

