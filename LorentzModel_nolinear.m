function Lorentz=LorentzModel_nolinear(x, freq)
% CJE LorentzModel
% Lorentzian  = (1/pi) * (hwhm) / (deltaf^2 + hwhm^2) (Wolfram)
% Peak height of Lorentzian = 4 / (pi*hwhm)
% This defnition of the Lorentzian has Area = 1 

area = x(1); 
hwhm = x(2);
f0 = x(3);
phase = x(4);
baseline0 = x(5);
baseline1 = 0; %no linear baseline (for Cr frame fit)

Absorption = 1/(2*pi) * area * ones(size(freq)) .* hwhm ./ ((freq-f0).^2 + hwhm.^2);
Dispersion = 1/(2*pi) * area * (freq-f0) ./ ((freq-f0).^2 + hwhm.^2);

Lorentz = Absorption*cos(phase) + Dispersion * sin(phase) ...
	  + baseline0 + baseline1*(freq-f0);

% from Lorentz component of LorentzGauss model
%F = (x(1)*ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1))  ... 
%    .* (exp(x(6)*(freq-x(3)).*(freq-x(3)))) ... % gaussian
%    + x(4)*(freq-x(3)) ... % linear baseline
%    +x(5); % constant baseline
