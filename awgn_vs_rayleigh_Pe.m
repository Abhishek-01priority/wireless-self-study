%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =================================================================
% Simulate Probability of Error for AWGN vs Rayleigh fading channel 
% for BPSK modulated signal
% =================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 1024;
% transmitted BPSK signal
signalTX = rand(1, N) - 0.5;
signalTX(find(signalTX <= 0)) = -1;
signalTX(find(signalTX > 0))  = +1;
sigPower = sum((signalTX.^2)/length(signalTX));

% snr
snr_in_dB = linspace(0.001,15,10);

for s = 1:length(snr_in_dB)
  
  % AWGN noise addition
  noisePwr = 10 ^ (-snr_in_dB(s)/10) * sigPower;
  noise = sqrt(noisePwr/2) .* randn(1,N);

  received_signal = signalTX + noise;

  % signal detection
  signalRX = received_signal;
  signalRX(find(signalRX <= 0)) = -1;
  signalRX(find(signalRX > 0)) = +1;
  %signalRX = floor(signalRX + 0.5);

  % Probability of Error calculation
  Pe_AWGN(s) = length(find(signalRX ~= signalTX)) / N;

end % for end

for s = 1:length(snr_in_dB)
  
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

% theoritical formula
Pe_theoritical_Fading = 0.5 * (1 - sqrt(snr_in_dB ./ (1 + snr_in_dB)));

% plot
semilogy(snr_in_dB, Pe_AWGN, '-ko', 'MarkerFaceColor', 'k')
hold on;
semilogy(snr_in_dB, Pe_Fading, '-r^', 'MarkerFaceColor', 'r')
title("Probability of Error in AWGN vs Rayleigh fading channel")
legend("awgn", "Rayleigh fading channel")
xlabel("SNR in dB")
ylabel("Probability of Error")
grid minor
