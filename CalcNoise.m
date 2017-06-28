function [MRS_struct, noise] = CalcNoise(MRS_struct, ii, subspec)
% Estimate noise in downfield frequency-domain signal

freq = MRS_struct.spec.freq;
if strcmp(subspec,'DIFF')
    spec = MRS_struct.spec.diff;
elseif strcmp(subspec,'OFF')
    spec = MRS_struct.spec.off;
elseif strcmp(subspec,'Water')
    spec = MRS_struct.spec.water;
end

N = 2; % fit second-order polynomial

indA = freq >= 10 & freq <= 11;
noiseA = real(spec(ii,indA));
pA = polyfit(freq(indA), noiseA, N);
noiseA_fit = polyval(pA, freq(indA));
noiseA_detrended = noiseA - noiseA_fit;

indB = freq >= 11 & freq <= 12;
noiseB = real(spec(ii,indB));
pB = polyfit(freq(indB), noiseB, N);
noiseB_fit = polyval(pB, freq(indB));
noiseB_detrended = noiseB - noiseB_fit;

sigma = min([std(noiseA_detrended) std(noiseB_detrended)]);
noise = 2*sigma;

eval(['MRS_struct.out.noise_' subspec '_A(ii,:) = noiseA;']);
eval(['MRS_struct.out.noise_' subspec '_A_detrended(ii,:) = noiseA_detrended;']);
eval(['MRS_struct.out.noise_' subspec '_B(ii,:) = noiseB;']);
eval(['MRS_struct.out.noise_' subspec '_B_detrended(ii,:) = noiseB_detrended;']);

% figure(332);
% clf;
% subaxis(2,1,1);
% hold on;
% plot(freq(indA), noiseA, 'r', freq(indA), noiseFitA, '--r');
% plot(freq(indB), noiseB, 'r', freq(indB), noiseFitB, '--r');
% hold off;
%
% minYlim = min([noiseA(:); noiseB(:)]);
% minYlim = minYlim - 0.25*abs(minYlim);
% maxYlim = max([noiseA(:); noiseB(:)]);
% maxYlim = maxYlim + 0.25*maxYlim;
%
% ylim([minYlim, maxYlim]);
% set(gca, 'xdir', 'reverse');
%
% subaxis(2,1,2);
% hold on;
% plot(freq(indA), noise_detrendedA, 'r');
% plot(freq(indB), noise_detrendedB, 'r');
% hold off;
%
% minYlim = min([min(noise_detrendedA) min(noise_detrendedB)]);
% minYlim = minYlim - 3*abs(minYlim);
% maxYlim = max([max(noise_detrendedA) max(noise_detrendedB)]);
% maxYlim = maxYlim + 3*maxYlim;
%
% ylim([minYlim, maxYlim]);
% set(gca, 'xdir', 'reverse');

