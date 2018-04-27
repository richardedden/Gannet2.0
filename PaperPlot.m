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
%   varargin: Optional inputs (entered as parameter-value pairs).
%                   nSpec: Number of spectra to plot, entered as a scalar.
%                          All spectra are plotted by default.
%                   freqLim: Limits of ppm axis, entered as a vector.
%                            Default is [0.5 4.5].
%                   signalLim: Limits of signal axis, entered as a vector.
%                              Default is [-6e-5 8e-5].
%                   signalLim2: Limits of signal axis for second edited
%                               metabolite if HERMES data.
%                   plotModel: Plot signal model fit, entered as a logical.
%                              Default is false.
%                   plotAvg: Plot averaged spectrum +/- 1 stdev, entered as
%                            a logical. Default is false.
%
% Examples:
%   PaperPlot(MRS_struct, 'nSpec', [1 3 4]);
%       This will plot the 1st, 3rd and 4th difference spectra in
%       MRS_struct along with the model fits of the peak(s) specified in
%       MRS_struct.p.target.
%
%   PaperPlot(MRS_struct, 'freqLim', [2.5 3.5], 'plotModel', true);
%       This will plot all difference spectra in MRS_struct with the model
%       fits of the peak(s) and limit the ppm axis from 2.5 to 3.5 ppm

% MM (180427)

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

% Check if data is PRIAM
if isfield(MRS_struct.p,'PRIAM')
    PRIAM = MRS_struct.p.PRIAM;
else
    PRIAM = 0;
end

% Set some defaults
vox = MRS_struct.p.Vox;
if ~PRIAM
    vox = vox(1);
end
defaultTarget = MRS_struct.p.target;
if HERMES
    defaultTarget2 = MRS_struct.p.target2;
    defaultnSpec = 1:size(MRS_struct.spec.(vox{1}).(defaultTarget).diff,1);
else
    defaultnSpec = 1:size(MRS_struct.spec.(vox{1}).(defaultTarget).diff,1);
end
defaultFreqLim = [0.5 4.5];
defaultSignalLim = [];
defaultSignalLim2 = [];
defaultPlotModel = false;
defaultPlotAvg = false;
expectedTargets = {'GABA', 'Glx', 'GABAGlx', 'GSH', 'Lac'};
grey = [0.6 0.6 0.6];
shading = 0.3;

% Parse input arguments
p = inputParser;
p.CaseSensitive = false;
p.addParamValue('target', defaultTarget, @(x) any(validatestring(x,expectedTargets))); %#ok<*NVREPL>
if HERMES
    p.addParamValue('target2', defaultTarget2, @(x) any(validatestring(x,expectedTargets)));
end
p.addParamValue('nSpec', defaultnSpec);
p.addParamValue('freqLim', defaultFreqLim);
p.addParamValue('signalLim', defaultSignalLim);
p.addParamValue('signalLim2', defaultSignalLim2);
p.addParamValue('plotModel', defaultPlotModel, @(x) islogical(x));
p.addParamValue('plotAvg', defaultPlotAvg, @(x) islogical(x));
p.parse(varargin{:});

target{1} = p.Results.target;
if HERMES
    target{2} = p.Results.target2;
end
nSpec = p.Results.nSpec;
freqLim = p.Results.freqLim;
signalLim = p.Results.signalLim;
signalLim2 = p.Results.signalLim2;
plotModel = p.Results.plotModel;
plotAvg = p.Results.plotAvg;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Plot spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq = MRS_struct.spec.freq;

if isfield(MRS_struct.spec.(vox{1}).(target{1}),'diff_scaled')
    diff = 'diff_scaled';
    fprintf('\nNB: Spectra are normalized to the amplitude of the respective modeled unsuppressed water reference signal.\n\n');
else
    diff = 'diff';
end

for ii = 1:length(vox)
    
    H = figure(199+ii);
    scr_sz = get(0, 'ScreenSize');
    fig_w = 1000;
    if HERMES
        fig_h = 1200;
    else
        fig_h = 707;
    end
    set(H, 'Color', 'w', 'Position', [(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
    clf;
    
    if HERMES
        subplot(2,1,1);
    end
    
    if isfield(MRS_struct.out.(vox{ii}),'water')
        scaleFactor = MRS_struct.out.(vox{ii}).water.ModelParam(:,1);
    else
        scaleFactor = ones(1,length(nSpec));
    end
    
    switch target{1}
        case 'GABA'
            modelFreq = freq(freq <= 3.55 & freq >= 2.79);
            
            % If nSpec > 1, find mean + stdevs
            if numel(nSpec) > 1 && plotAvg
                mu = mean(MRS_struct.spec.(vox{ii}).GABA.(diff),1);
                sigma = std(MRS_struct.spec.(vox{ii}).GABA.(diff),0,1);
                UB = mu + sigma;
                LB = mu - sigma;
            end
            
            if numel(nSpec) > 1 && plotAvg
                hold on;
                h = patch([freq fliplr(freq)], [real(UB) fliplr(real(LB))], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
                plot(freq, real(mu), 'k');
                hold off;
            else
                hold on;
                for jj = 1:length(nSpec)
                    if plotModel
                        h(:,jj) = plot(freq, real(MRS_struct.spec.GABAGlx.(diff)(nSpec(jj),:)), 'k', ...
                            modelFreq, GABAGlxModel(MRS_struct.out.GABA.ModelFit(nSpec(jj),:),modelFreq) ./ scaleFactor(nSpec(jj)), 'r', 'LineWidth', 1);
                    else
                        h = plot(freq, real(MRS_struct.spec.GABAGlx.(diff)(nSpec(jj),:)), 'k', 'LineWidth', 1);
                    end
                end
                hold off;
            end
            
        case 'GABAGlx'
            modelFreq = freq(freq <= 4.1 & freq >= 2.79);
            
            % If nSpec > 1, find mean + stdevs
            if numel(nSpec) > 1 && plotAvg
                mu = mean(MRS_struct.spec.(vox{ii}).GABAGlx.(diff),1);
                sigma = std(MRS_struct.spec.(vox{ii}).GABAGlx.(diff),0,1);
                UB = mu + sigma;
                LB = mu - sigma;
            end
            
            if numel(nSpec) > 1 && plotAvg
                hold on;
                h = patch([freq fliplr(freq)], [real(UB) fliplr(real(LB))], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
                plot(freq, real(mu), 'k');
                hold off;
            else
                hold on;
                for jj = 1:length(nSpec)
                    if plotModel
                        h(:,jj) = plot(freq, real(MRS_struct.spec.(vox{ii}).GABAGlx.(diff)(nSpec(jj),:)), 'k', ...
                            modelFreq, GABAGlxModel(MRS_struct.out.(vox{ii}).GABA.ModelParam(nSpec(jj),:),modelFreq) ./ scaleFactor(nSpec(jj)), 'r', 'LineWidth', 1);
                    else
                        h(jj) = plot(freq, real(MRS_struct.spec.(vox{ii}).GABAGlx.(diff)(nSpec(jj),:)), 'k', 'LineWidth', 1);
                    end
                end
                hold off;
            end
            
        case 'GSH'
            modelFreq = freq(freq <= 3.3 & freq >= 2.25);
            GSHgaussModel = eval(['@' MRS_struct.p.GSH_model 'Model']);
            
            % If nSpec > 1, find mean + stdevs
            if numel(nSpec) > 1 && plotAvg
                mu = mean(MRS_struct.spec.(vox{ii}).GSH.(diff),1);
                sigma = std(MRS_struct.spec.(vox{ii}).GSH.(diff),0,1);
                UB = mu + sigma;
                LB = mu - sigma;
            end
            
            if numel(nSpec) > 1 && plotAvg
                hold on;
                patch([freq fliplr(freq)], [real(UB) fliplr(real(LB))], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
                h = plot(freq, real(mu), 'k');
                hold off;
            else
                hold on;
                for jj = 1:length(nSpec)
                    if plotModel
                        h(:,jj) = plot(freq, real(MRS_struct.spec.(vox{ii}).GSH.(diff)(nSpec(jj),:)), 'k', ...
                            modelFreq, GSHgaussModel(MRS_struct.out.(vox{ii}).GSH.ModelParam(nSpec(jj),:),modelFreq) ./ scaleFactor(nSpec(jj)), 'r', 'LineWidth', 1);
                    else
                        h = plot(freq, real(MRS_struct.spec.(vox{ii}).GSH.(diff)(nSpec(jj),:)), 'k', 'LineWidth', 1);
                    end
                end
                hold off;
            end
    end
    
    set(gca, 'TickDir', 'out', 'XLim', freqLim, 'XDir', 'reverse', 'YTick', [], 'YColor', 'w', 'Box', 'off', 'FontSize', 20, 'LineWidth',1);
    xlabel('ppm', 'FontWeight', 'bold', 'FontSize', 28);
    
    % Set YLim
    if isempty(signalLim)
        switch target{1}
            case {'GABA','GABAGlx','Glx'}
                peakrange = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.24;
            case 'GSH'
                peakrange = MRS_struct.spec.freq <= 3.5 & MRS_struct.spec.freq >= 1;
            case 'Lac'
                peakrange = MRS_struct.spec.freq <= 4.5 & MRS_struct.spec.freq >= 0.5;
        end
        for jj = 1:size(h,2)
            maxPeakHeight(jj) = abs(max(h(1,jj).YData(peakrange))); %#ok<*AGROW>
        end
        yaxismax = max(maxPeakHeight);
        yaxismax = yaxismax + yaxismax/10;
        yaxismin = -yaxismax/2;
        
        if isfield(MRS_struct.spec.(vox{1}).(target{1}),'diff_scaled')
            signalLim = [-6e-5 8e-5];
        else
            signalLim = [yaxismin yaxismax];
        end
    end
    set(gca, 'YLim', signalLim);
    
    % For HERMES data
    if HERMES
        subplot(2,1,2);
        switch target{2}
            case 'GSH'
                modelFreq = freq(freq <= 3.3 & freq >= 2.25);
                GSHgaussModel = eval(['@' MRS_struct.p.GSH_model 'Model']);
                
                % If nSpec > 1, find mean + stdevs
                if numel(nSpec) > 1 && plotAvg
                    mu = mean(MRS_struct.spec.(vox{ii}).GSH.(diff),1);
                    sigma = std(MRS_struct.spec.(vox{ii}).GSH.(diff),0,1);
                    UB = mu + sigma;
                    LB = mu - sigma;
                end
                
                if numel(nSpec) > 1 && plotAvg
                    hold on;
                    h2 = patch([freq fliplr(freq)], [real(UB) fliplr(real(LB))], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
                    plot(freq, real(mu), 'k');
                    hold off;
                else
                    hold on;
                    for jj = 1:length(nSpec)
                        if plotModel
                            h2(:,jj) = plot(freq, real(MRS_struct.spec.(vox{ii}).GSH.(diff)(nSpec(jj),:)), 'k', ...
                                modelFreq, GSHgaussModel(MRS_struct.out.(vox{ii}).GSH.ModelParam(nSpec(jj),:),modelFreq) ./ scaleFactor(nSpec(jj)), 'r', 'LineWidth', 1);
                        else
                            h2 = plot(freq, real(MRS_struct.spec.(vox{ii}).GSH.(diff)(nSpec(jj),:)), 'k', 'LineWidth', 1);
                        end
                    end
                    hold off;
                end
                
%             case 'Lac'
%                 lb = find(specfreq <= 4.10);
%                 ub = find(specfreq >= 2.79);
%                 range = intersect(lb,ub);
%                 modelfreq = specfreq(range);
%
%                 for ii = 1:length(p.Results.nSpec)
%                     hold on;
%                     plot(specfreq, real(MRS_struct.spec.(diff)(p.Results.nSpec(ii),:)), 'k', ...
%                         modelfreq, GABAGlxModel(MRS_struct.out.GABA.ModelFit(p.Results.nSpec(ii),:),modelfreq), 'r');
%                     hold off;
%                 end
        end
        
        set(gca, 'TickDir', 'out', 'XLim', freqLim, 'XDir', 'reverse', 'YLim', signalLim, 'YTick', [], 'YColor', 'w', 'Box', 'off', 'FontSize', 20, 'LineWidth',1);
        xlabel('ppm', 'FontWeight', 'bold', 'FontSize', 28);
        
        % Set YLim
        if isempty(signalLim2)
            switch target{2}
                case {'GABA','GABAGlx','Glx'}
                    peakrange = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.24;
                case 'GSH'
                    peakrange = MRS_struct.spec.freq <= 3.5 & MRS_struct.spec.freq >= 1;
                case 'Lac'
                    peakrange = MRS_struct.spec.freq <= 4.5 & MRS_struct.spec.freq >= 0.5;
            end
            for jj = 1:size(h2,2)
                maxPeakHeight(jj) = abs(max(h2(1,jj).YData(peakrange)));
            end
            yaxismax = max(maxPeakHeight);
            yaxismax = yaxismax + yaxismax/4;
            yaxismin = -yaxismax*2;
            
            if isfield(MRS_struct.spec.(vox{1}).(target{2}),'diff_scaled')
                signalLim2 = [-10e-5 4e-5];
            else
                signalLim2 = [yaxismin yaxismax];
            end
        end
        set(gca, 'YLim', signalLim2);
        
    end
    
end




