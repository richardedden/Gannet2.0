function PaperPlot(MRS_struct, varargin)
% PaperPlot(MRS_struct, varargin)
%
% This function will plot the difference spectrum, or any number of
% difference spectra, saved in MRS_struct along with the corresponding
% model fits of GABA/Glx, GSH and/or Lac peak(s) (as specified by
% MRS_struct.p.target and .target2). You can choose to plot a single
% spectrum, a select number of spectra or all spectra. Multiple spectra
% will be plotted in the same figure. If data were acquired with HERMES,
% then the relevant edited spectra will be plotted in separate subplots.
%
% Inputs:
%   MRS_struct: Structure output from GannetFit (required).
%   varargin:   Optional inputs (entered as parameter-value pairs).
%                   nSpec: Number of spectra to plot, entered as a scalar.
%                          All spectra are plotted by default.
%                   freqLim: Limits of ppm axis, entered as a vector.
%                            Default is [0.5 4.5].
%                   plotAvg: Plot averaged spectrum +/- 1 stdev, entered as
%                            a logical. Default is false.
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

% MM (161024)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Parse inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    error('Not enough inputs! MRS_struct is required!');
end

% Check if data is MEGA- or HERMES-edited
if isfield(MRS_struct.p,'HERMES')
    HERMES = MRS_struct.p.HERMES;
else
    HERMES = 0;
end

% Set some defaults
defaultTarget = MRS_struct.p.target;
if HERMES; defaultTarget2 = MRS_struct.p.target2; end
if HERMES
    defaultnSpec = eval(['1:size(MRS_struct.spec.' defaultTarget '.diff,1)']);
else
    defaultnSpec = 1:size(MRS_struct.spec.diff,1);
end
defaultFreqLim = [0.5 4.5];
defaultSignalLim = [-0.05 0.03];
defaultPlotAvg = false;
expectedTargets = {'GABA', 'Glx', 'GABAGlx', 'GSH'};
grey = [0.6 0.6 0.6];
shading = 0.3;

% Parse input arguments
p = inputParser;
p.CaseSensitive = false;
p.addParamValue('target', defaultTarget, @(x) any(validatestring(x,expectedTargets))); %#ok<*NVREPL>
if HERMES; p.addParamValue('target2', defaultTarget2, @(x) any(validatestring(x,expectedTargets))); end
p.addParamValue('nSpec', defaultnSpec);
p.addParamValue('freqLim', defaultFreqLim);
p.addParamValue('signalLim', defaultSignalLim);
p.addParamValue('plotAvg', defaultPlotAvg, @(x) islogical(x));
p.parse(varargin{:});

target = p.Results.target;
if HERMES; target2 = p.Results.target2; end
nSpec = p.Results.nSpec;
freqLim = p.Results.freqLim;
signalLim = p.Results.signalLim;
plotAvg = p.Results.plotAvg;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Plot spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

specFreq = MRS_struct.spec.freq(1,:);
h = figure(200);
set(h, 'Color', 'w', 'Units', 'Normalized', 'OuterPosition', [0 0.4 1 1-0.4]);
clf;

if HERMES
    subplot(1,2,1);
end

switch target
    case 'GABA'
        lb = find(specFreq <= 3.55);
        ub = find(specFreq >= 2.79);
        range = intersect(lb,ub);
        modelFreq = specFreq(range);
        
        % If nSpec > 1, find mean + stdevs
        if numel(nSpec) > 1 && plotAvg
            m = mean(MRS_struct.spec.diff,1);
            s = std(MRS_struct.spec.diff,0,1);
            sUB = m + s;
            sLB = m - s;
        end
        
        if plotAvg
            hold on;
            patch([specFreq fliplr(specFreq)], [sUB fliplr(sLB)], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
            plot(specFreq, m, 'k');
            hold off;
        else
            hold on;
            for ii = nSpec
                plot(specFreq, MRS_struct.spec.GABAGlx.diff(nSpec(ii),:), 'k', ...
                    modelFreq, GABAGlxModel(MRS_struct.out.GABA.ModelFit(nSpec(ii),:),modelFreq), 'r');
            end
            hold off;
        end
        
    case 'GABAGlx'
        lb = find(specFreq <= 4.10);
        ub = find(specFreq >= 2.79);
        range = intersect(lb,ub);
        modelFreq = specFreq(range);
        
        % If nSpec > 1, find mean + stdevs
        if numel(nSpec) > 1 && plotAvg
            m = mean(MRS_struct.spec.GABAGlx.diff,1);
            s = std(MRS_struct.spec.GABAGlx.diff,0,1);
            sUB = m + s;
            sLB = m - s;
        end
        
        if plotAvg
            hold on;
            patch([specFreq fliplr(specFreq)], [sUB fliplr(sLB)], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
            plot(specFreq, m, 'k');
            hold off;
        else
            for ii = nSpec
                hold on;
                plot(specFreq, MRS_struct.spec.GABAGlx.diff(nSpec(ii),:), 'k', ...
                    modelFreq, GABAGlxModel(MRS_struct.out.GABA.ModelFit(nSpec(ii),:),modelFreq), 'r');
                hold off;
            end
        end
        
    case 'GSH'
%         lb = find(specfreq <= 4.10);
%         ub = find(specfreq >= 2.79);
%         range = intersect(lb,ub);
%         modelfreq = specfreq(range);
%         
%         for ii = 1:length(p.Results.nSpec)
%             hold on;
%             plot(specfreq, MRS_struct.spec.diff(p.Results.nSpec(ii),:), 'k', ...
%                 modelfreq, GABAGlxModel(MRS_struct.out.GABA.ModelFit(p.Results.nSpec(ii),:),modelfreq), 'r');
%             hold off;
%         end
end

set(gca, 'TickDir', 'out', 'XLim', freqLim, 'XDir', 'reverse', 'YLim', signalLim, 'YTick', [], 'YColor', 'w', 'Box', 'off');
xlabel('ppm', 'FontWeight', 'bold');

% For HERMES data
if HERMES    
    subplot(1,2,2);
    switch target2
        case 'GSH'
%             lb = find(specfreq <= 4.10);
%             ub = find(specfreq >= 2.79);
%             range = intersect(lb,ub);
%             modelfreq = specfreq(range);

        % If nSpec > 1, find mean + stdevs
        if numel(nSpec) > 1 && plotAvg
            m = mean(MRS_struct.spec.GSH.diff,1);
            s = std(MRS_struct.spec.GSH.diff,0,1);
            sUB = m + s;
            sLB = m - s;
        end
        
        if plotAvg
            hold on;
            patch([specFreq fliplr(specFreq)], [sUB fliplr(sLB)], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
            plot(specFreq, m, 'k');
            hold off;
        else
            for ii = 1:length(nSpec)
                hold on;
                plot(specFreq, MRS_struct.spec.GSH.diff(nSpec(ii),:), 'k'); % MM: need to add model fit
                hold off;
            end
        end
            
            freqLim = [0.5 4.0];
            
        case 'Lac'
%             lb = find(specfreq <= 4.10);
%             ub = find(specfreq >= 2.79);
%             range = intersect(lb,ub);
%             modelfreq = specfreq(range);
%             
%             for ii = 1:length(p.Results.nSpec)
%                 hold on;
%                 plot(specfreq, MRS_struct.spec.diff(p.Results.nSpec(ii),:), 'k', ...
%                     modelfreq, GABAGlxModel(MRS_struct.out.GABA.ModelFit(p.Results.nSpec(ii),:),modelfreq), 'r');
%                 hold off;
%             end
    end
    
    set(gca, 'TickDir', 'out', 'XLim', freqLim, 'XDir', 'reverse', 'YLim', signalLim.*[0.5 2.5], 'YTick', [], 'YColor', 'w', 'Box', 'off');
    xlabel('ppm', 'FontWeight', 'bold');
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Print spectra (optional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







