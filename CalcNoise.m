function noise = CalcNoise(freq, spec)
% Estimate noise in downfield frequency-domain signal

N = 2; % fit second-order polynomial

indA = freq >= 10 & freq <= 11;
noiseA = real(spec(indA));
pA = polyfit(freq(indA), noiseA, N);
noiseA_fit = polyval(pA, freq(indA));
noiseA_detrended = noiseA - noiseA_fit;

indB = freq >= 11 & freq <= 12;
noiseB = real(spec(indB));
pB = polyfit(freq(indB), noiseB, N);
noiseB_fit = polyval(pB, freq(indB));
noiseB_detrended = noiseB - noiseB_fit;

sigma = min([std(noiseA_detrended) std(noiseB_detrended)]);
noise = 2*sigma;
