function PaperPlot(MRS_struct, varargin)
% PaperPlot(MRS_struct, varargin)
%
% This function will plot the difference spectrum, or any number of
% difference spectra, saved in MRS_struct along with the corresponding
% model fits of the GABA, Glx or GABA/Glx peak(s) (as specified by
% MRS_struct.p.target). You can choose to plot a single spectrum, a select
% number of spectra or all spectra. Multiple spectra will be plotted in
% the same figure.
%
% For paper output, save in .eps format (matlab .pdf isn't good).
%
% Inputs:
%   MRS_struct: Structure output from GannetFit (required).
%   varargin:   Optional inputs (entered as parameter-value pairs).
%                   nSpec: Number of spectra to plot, entered as a vector.
%                          All spectra are plotted by default.
%                   freqLim: Limits of ppm axis. Default is [0.5 4.5].
%
% Examples:
%   PaperPlot(MRS_struct, 'nSpec', [1 3 4];
%       This will plot the 1st, 3rd and 4th difference spectra in
%       MRS_struct along with the model fits of the peak(s) specified in
%       MRS_struct.p.target.
%
%   PaperPlot(MRS_struct, 'freqLim', [2.5 3.5]);
%       This will plot all difference spectra in MRS_struct along with the
%       model fits of the peak(s) specified in MRS_struct.p.target and
%       limit the ppm axis from 2.5 to 3.5 ppm

% MM (160720)

%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Parse inputs
%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    error('Too few inputs! MRS_struct is required!');
end

% Set some defaults
defaultTarget = MRS_struct.p.target;
defaultnSpec = 1:size(MRS_struct.spec.diff,1);
defaultFreqLim = [0.5 4.5];
expectedTargets = {'GABA', 'Glx', 'GABAGlx'};

% Parse argmuments
p = inputParser;
p.addParamValue('target', defaultTarget, @(x) any(validatestring(x,expectedTargets)));
p.addParamValue('nSpec', defaultnSpec);
p.addParamValue('freqLim', defaultFreqLim);
p.parse(varargin{:});


%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Plot spectra
%%%%%%%%%%%%%%%%%%%%%%%%

specfreq = MRS_struct.spec.freq(1,:);
figure('Color', 'w', 'Units', 'Normalized', 'OuterPosition', [0 0.05 0.65 1-0.05]);

switch p.Results.target
    case 'GABA'
        lb = find(specfreq <= 3.55);
        ub = find(specfreq >= 2.79);
        range = intersect(lb,ub);
        modelfreq = specfreq(range);
        
        for ii = 1:length(p.Results.nSpec)
            hold on;
            plot(specfreq, MRS_struct.spec.diff(p.Results.nSpec(ii),:), 'k', ...
                modelfreq, GaussModel(MRS_struct.out.GABA.ModelFit(p.Results.nSpec(ii),:),modelfreq), 'r');
            hold off;
        end
                
    case 'Glx'
        lb = find(specfreq <= 4.10);
        ub = find(specfreq >= 3.45);
        range = intersect(lb,ub);
        modelfreq = specfreq(range);
        
        for ii = 1:length(p.Results.nSpec)
            hold on;
            plot(specfreq, MRS_struct.spec.diff(p.Results.nSpec(ii),:), 'k', ...
                modelfreq, DoubleGaussModel(MRS_struct.out.Glx.ModelFit(p.Results.nSpec(ii),:),modelfreq), 'r');
            hold off;
        end
                
    case 'GABAGlx' 
        lb = find(specfreq <= 4.10);
        ub = find(specfreq >= 2.79);
        range = intersect(lb,ub);
        modelfreq = specfreq(range);
        
        for ii = 1:length(p.Results.nSpec)
            hold on;
            plot(specfreq, MRS_struct.spec.diff(p.Results.nSpec(ii),:), 'k', ...
                modelfreq, GABAGlxModel(MRS_struct.out.GABA.ModelFit(p.Results.nSpec(ii),:),modelfreq), 'r');
            hold off;
        end        
end

set(gca, 'XLim', p.Results.freqLim, 'XDir', 'reverse', 'YTick', [], 'YColor', 'w', 'Box', 'off');
xlabel('ppm', 'FontWeight', 'bold');





