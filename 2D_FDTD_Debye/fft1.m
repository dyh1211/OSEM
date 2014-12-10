function [f, Y] = fft1(y, Fs)
	L = length(y);
	NFFT = 2^nextpow2(L); % Next power of 2 from length of y
	Y = fft(y, NFFT) / L;
	f = Fs / 2 * linspace(0, 1, NFFT / 2 + 1);
	Y = 2 * abs(Y(1:NFFT / 2 + 1));
	Y = Y ./ max(Y);
end

