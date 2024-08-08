function [Y] = dbaf(x, N)
	Y = 10*log10( abs ( fftshift ( fft (x, N) ) ) );
end
