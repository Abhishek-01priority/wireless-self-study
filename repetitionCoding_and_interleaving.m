%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =================================================================
% This file investigates the effect of interleaving and repetition 
% coding on a fading channel.
% TODO: Move AWGN and fading randn outside with fixed seed for fair 
% comparision
% =================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixing the randn seed.
V = randn("seed");

% generate clarkes model
% 22 samples in doppler spectrum and 100Hz of doppler shift
[h, number_of_samples] = clarke_model(22, 100); 

N = number_of_samples;
% transmitted BPSK signal
signalTX = rand(1, N) - 0.5;
signalTX(find(signalTX <= 0)) = -1;
signalTX(find(signalTX > 0))  = +1;
sigPower = sum((signalTX.^2)/length(signalTX));

% snr
snr_in_dB = linspace(-15,15,10);

% repetition coding: 
% Each sample is repeated reptitionCodingRate times.
% So total size of transmitted signal is reptitionCodingRate*N
reptitionCodingRate = 6;
Repeated_SignalTX = reshape( repmat(signalTX, reptitionCodingRate, 1), 1, [] );

% interleaving after repetition coding
% the repeated symbols are spread accross coherence time.
Interleaved_SignalTX = repmat(signalTX, 1, reptitionCodingRate);

N_after_rep = N * reptitionCodingRate;

% channel simulation for repetition and interleaving
for s = 1:length(snr_in_dB)
  
  randn("seed", V);
  % AWGN noise addition. transmited signal power for repetition and interleaving is
  % equal to transmited power of non-repeated-non-interleaved signal.
  noisePwr = 10 ^ (-snr_in_dB(s)/10) * sigPower;
  noise = sqrt(noisePwr/2) .* ((randn(1,N_after_rep)) + 1i * (randn(1,N_after_rep)));
  noise = noise - mean(noise);

  % Fading channel:
  % As such, the fading channel is modelled as complex gaussian with 0 mean and 1 variance
  % norm(h) is rayleigh distributed and square(norm(h)) is exponential
  
  % TODO: understand when clarkes model is used. 
  %h = clarke_model(N, 100, 2);
  %received_signal = ifft (fft(h,N) .* fft(signalTX, N), N) + noise;

  % h = CN(0, 1)
  %h = sqrt(1/2) * ( randn(1,N_after_rep) + 1i * randn(1,N_after_rep));
  %h = h - mean(h);

  % Idea is to ensure that repetition coding and interleaving undergo same 
  % AWGN noise and fading
  rx_repeated_signal = h .* Repeated_SignalTX  + noise;
  rx_interleaved_signal = h .* Interleaved_SignalTX + noise;

  % signal estimation
  y_tilda_repetition = conj(h) .* rx_repeated_signal ./ sqrt(conj(h) .* h);
  y_tilda_interleaved = conj(h) .* rx_interleaved_signal ./ sqrt(conj(h) .* h);

  % coherent accumulation
  coherent_y_rep = sum ( reshape( y_tilda_repetition, reptitionCodingRate, [] ) ); % row-wise sum
  coherent_y_int = sum ( reshape( y_tilda_interleaved, N, [] ) , 2)'; % columnwise sum

  % signal detection
  signalRX_repetition = real(coherent_y_rep);
  signalRX_repetition(find(signalRX_repetition <= 0)) = -1;
  signalRX_repetition(find(signalRX_repetition > 0)) = +1;

  signalRX_interleaved = real(coherent_y_int);
  signalRX_interleaved(find(signalRX_interleaved <= 0)) = -1;
  signalRX_interleaved(find(signalRX_interleaved > 0)) = +1;

  % Probability of Error calculation
  Pe_repetition(s) = length(find(signalRX_repetition ~= signalTX)) / N;
  Pe_interleaved(s) = length(find(signalRX_interleaved ~= signalTX)) / N;

end % for end

% channel sim for non-repetition and non-interleaving
for s = 1:length(snr_in_dB)
  
  randn("seed", V);
  % AWGN noise addition
  noisePwr = 10 ^ (-snr_in_dB(s)/10) * sigPower;
  noise = sqrt(noisePwr/2) .* ((randn(1,N)) + 1i * (randn(1,N)));

  % Fading channel:
  % TODO: understand when clarkes model is used. 
  % As such, the fading channel is modelled as complex gaussian with 0 mean and 1 variance
  % norm(h) is rayleigh distributed and square(norm(h)) is exponential
  
  %h = clarke_model(N, 100, 2);
  %received_signal = ifft (fft(h,N) .* fft(signalTX, N), N) + noise;
  h = sqrt(1/2) * ( randn(1,N) + 1i * randn(1,N));
  h = h - mean(h);
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
%semilogy(snr_in_dB, Pe_Fading, '-bs', 'MarkerFaceColor', 'b')
%title("Probability of Error in AWGN vs Rayleigh fading channel")
legend("Repetition coding", "Repetition and Interleaving", "No repetition and no interleaving")
xlabel("SNR in dB")
ylabel("Probability of Error")
grid minor


