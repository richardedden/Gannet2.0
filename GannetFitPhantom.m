function MRS_struct = GannetFitPhantom(MRS_struct, varargin)
% Gannet 3.0 GannetFitPhantom
% Updates by MM 2018

if MRS_struct.p.PRIAM
    vox = MRS_struct.p.Vox;
else
    vox = {MRS_struct.p.Vox{1}};
end

if MRS_struct.p.HERMES && nargin < 2
    target = {MRS_struct.p.target, MRS_struct.p.target2};
else
    % varargin = Optional arguments if user wants to overwrite fitting
    %            parameters set in GannetPreInitialise; can include several
    %            options, which are:
    %            'GABA' or 'Glx': target metabolite
    if nargin > 1
        switch varargin{1}
            case 'GABA'
                MRS_struct.p.target = 'GABA';
            case 'Glx'
                MRS_struct.p.target = 'Glx';
            case 'GSH'
                MRS_struct.p.target = 'GSH';
            case 'Lac'
                MRS_struct.p.target = 'Lac';
            case 'EtOH'
                MRS_struct.p.target = 'EtOH';
        end
    end
    target = {MRS_struct.p.target};
end

freq = MRS_struct.spec.freq;
MRS_struct.version.fit = '180912';

lsqopts = optimset('lsqcurvefit');
lsqopts = optimset(lsqopts,'MaxIter',1e5,'MaxFunEvals',1e5,'TolX',1e-10,'TolFun',1e-10,'Display','off');
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts,'MaxIter',1e5,'MaxFunEvals',1e5,'TolX',1e-10,'TolFun',1e-10);

% Loop over voxels if PRIAM
for kk = 1:length(vox)
    
    if strcmp(MRS_struct.p.Reference_compound,'H2O')
        WaterData = MRS_struct.spec.(vox{kk}).water;
    end
    
    % Loop over edited spectra if HERMES
    for trg = 1:length(target)
        
        fprintf('\nFitting %s...',target{trg});
        
        % Defining variables -- MGSaleh 2016
        DIFF = MRS_struct.spec.(vox{kk}).(target{trg}).diff;
        numscans = size(DIFF,1);
        
        for ii = 1:numscans
            
            if strcmp(target{trg},'GABA')
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   1.  GABA Fit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                freqbounds = find(freq <= 3.3 & freq >= 2.7);
                plotbounds = find(freq <= 3.4 & freq >= 2.6);
                
                maxinGABA = max(real(DIFF(ii,freqbounds)));
                grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                LinearInit = grad_points ./ abs(freq(1) - freq(2));
                T2 = 60;
                omega0 = 3.01;
                J = 7.5/MRS_struct.p.LarmorFreq(ii);
                
                ThreeLorentzModelInit = [maxinGABA/T2 maxinGABA/T2 maxinGABA/T2/3 ...
                                         T2 omega0 J 0 0 0 -LinearInit 0];
                lb = [-4000*maxinGABA/T2 -4000*maxinGABA/T2 -4000*maxinGABA/T2/3 ...
                      0 omega0-J J-0.01 -pi -pi -pi -40*maxinGABA -2000*maxinGABA];
                ub = [4000*maxinGABA/T2 4000*maxinGABA/T2 4000*maxinGABA/T2/3 ...
                      T2*100 omega0+J J+0.01 pi pi pi 40*maxinGABA 1000*maxinGABA];
                
                % Least-squares model fitting
                ThreeLorentzModelInit = lsqcurvefit(@ThreeLorentzModel, ThreeLorentzModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)), lb, ub, lsqopts);
                [ThreeLorentzModelParam, resid] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)), @ThreeLorentzModel, ThreeLorentzModelInit, nlinopts);
                
                GABAheight = max(ThreeLorentzModelParam(1:3)) * ThreeLorentzModelParam(4);
                MRS_struct.out.(vox{kk}).(target{trg}).FitError(ii) = 100*std(resid)/GABAheight;
                MRS_struct.out.(vox{kk}).(target{trg}).Area(ii) = sum(ThreeLorentzModel(ThreeLorentzModelParam,freq(freqbounds))) * abs(freq(1) - freq(2));
                MRS_struct.out.(vox{kk}).(target{trg}).FWHM(ii) = 1./(pi*ThreeLorentzModelParam(4))*1e3;
                MRS_struct.out.(vox{kk}).(target{trg}).ModelParam(ii,:) = ThreeLorentzModelParam;
                MRS_struct.out.(vox{kk}).(target{trg}).Resid(ii,:) = resid;
                
                % Calculate SNR of GABA signal (MM: 170502)
                noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                MRS_struct.out.(vox{kk}).(target{trg}).SNR(ii) = abs(GABAheight)/noiseSigma_DIFF;
                
            elseif strcmp(target{trg},'Glx')
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   1.  Glx Fit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                freqbounds = find(freq <= 4.05 & freq >= 3.45);
                plotbounds = find(freq <= 4.15 & freq >= 3.35);
                
                maxinGlx = max(real(DIFF(ii,freqbounds)));
                grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                LinearInit = grad_points ./ abs(freq(1) - freq(2));
                T2 = 60;
                omega0 = 3.74;
                J = 5.8/MRS_struct.p.LarmorFreq(ii);
                
                TwoLorentzModelInit = [maxinGlx/T2 maxinGlx/T2 ...
                                       T2 omega0 J 0 0 -LinearInit 0];
                lb = [-4000*maxinGlx/T2 -4000*maxinGlx/T2 ...
                      0 omega0-J J-0.01 -pi -pi -40*maxinGlx -2000*maxinGlx];
                ub = [4000*maxinGlx/T2 4000*maxinGlx/T2 ...
                      T2*100 omega0+J J+0.01 pi pi 40*maxinGlx 1000*maxinGlx];
                
                % Least-squares model fitting
                TwoLorentzModelInit = lsqcurvefit(@TwoLorentzModel2, TwoLorentzModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)), lb, ub, lsqopts);
                [TwoLorentzModelParam, resid] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)), @TwoLorentzModel2, TwoLorentzModelInit, nlinopts);
                
                Glxheight = max(TwoLorentzModelParam(1:2)) * TwoLorentzModelParam(3);
                MRS_struct.out.(vox{kk}).(target{trg}).FitError(ii) = 100*std(resid)/Glxheight;
                MRS_struct.out.(vox{kk}).(target{trg}).Area(ii) = sum(TwoLorentzModel2(TwoLorentzModelParam,freq(freqbounds))) * abs(freq(1) - freq(2));
                MRS_struct.out.(vox{kk}).(target{trg}).FWHM(ii) = 1./(pi*TwoLorentzModelParam(3))*1e3;
                MRS_struct.out.(vox{kk}).(target{trg}).ModelParam(ii,:) = TwoLorentzModelParam;
                MRS_struct.out.(vox{kk}).(target{trg}).Resid(ii,:) = resid;
                
                % Calculate SNR of Glx signal (MM: 170502)
                noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                MRS_struct.out.(vox{kk}).(target{trg}).SNR(ii) = abs(Glxheight)/noiseSigma_DIFF;
                
            elseif strcmp(target{trg},'GSH')
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   1.  GSH Fit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                freqbounds = find(freq <= 3.2 & freq >= 2.7);
                plotbounds = find(freq <= 3.4 & freq >= 2.5);
                
                maxinGSH = max(real(DIFF(ii,freqbounds)));
                grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                LinearInit = grad_points ./ abs(freq(1) - freq(2));
                T2 = 60;
                omega0 = 2.95;
                J = 4/MRS_struct.p.LarmorFreq(ii);
                
                SixLorentzModelInit = [-maxinGSH/T2/10 -maxinGSH/T2/10 ...
                                       maxinGSH/T2 maxinGSH/T2 ...
                                       -maxinGSH/T2/10 -maxinGSH/T2/10 ...
                                       T2 omega0 J ...
                                       180 180 0 0 180 180 ...
                                       -LinearInit 0];  
                                     
                lb = [-4000*maxinGSH/T2/10 -4000*maxinGSH/T2/10 ...
                      -4000*maxinGSH/T2    -4000*maxinGSH/T2 ...
                      -4000*maxinGSH/T2/10 -4000*maxinGSH/T2/10 ...
                      0 omega0-J J-0.01 ...
                      -pi -pi -pi -pi -pi -pi ...
                      -40*maxinGSH -2000*maxinGSH];
                  
                ub = [4000*maxinGSH/T2/10 4000*maxinGSH/T2/10 ...
                      4000*maxinGSH/T2    4000*maxinGSH/T2 ...
                      4000*maxinGSH/T2/10 4000*maxinGSH/T2/10 ...
                      T2*100 omega0+J J+0.01 ...
                      pi pi pi pi pi pi ...
                      40*maxinGSH 2000*maxinGSH];
                
                % Least-squares model fitting
                SixLorentzModelInit = lsqcurvefit(@SixLorentzModel, SixLorentzModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)), lb, ub, lsqopts);
                [SixLorentzModelParam, resid] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)), @SixLorentzModel, SixLorentzModelInit, nlinopts);
                
                GSHheight = max(SixLorentzModelParam(3:4)) * SixLorentzModelParam(7);
                MRS_struct.out.(vox{kk}).(target{trg}).FitError(ii) = 100*std(resid)/GSHheight;
                MRS_struct.out.(vox{kk}).(target{trg}).Area(ii) = sum(SixLorentzModel(SixLorentzModelParam,freq(freqbounds))) * abs(freq(1) - freq(2));
                MRS_struct.out.(vox{kk}).(target{trg}).FWHM(ii) = 1./(pi*SixLorentzModelParam(7))*1e3;
                MRS_struct.out.(vox{kk}).(target{trg}).ModelParam(ii,:) = SixLorentzModelParam;
                MRS_struct.out.(vox{kk}).(target{trg}).Resid(ii,:) = resid;
                
                % Calculate SNR of GABA signal (MM: 170502)
                noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                MRS_struct.out.(vox{kk}).(target{trg}).SNR(ii) = abs(GSHheight)/noiseSigma_DIFF;
                
            elseif strcmp(target{trg},'EtOH')
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   1.  EtOH Fit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                freqbounds = find(freq <= 1.4 & freq >= 0.9);
                plotbounds = find(freq <= 1.6 & freq >= 0.7);
                
                maxinEtOH = max(real(DIFF(ii,freqbounds)));
                grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                LinearInit = grad_points ./ abs(freq(1) - freq(2));
                T2 = 70;
                omega0 = 1.18;
                J = 7.5/MRS_struct.p.LarmorFreq(ii);
                
                ThreeLorentzModelInit = [maxinEtOH/T2 maxinEtOH/T2 maxinEtOH/T2/3 ...
                                         T2 omega0 J 0 0 0 -LinearInit 0];
                lb = [-4000*maxinEtOH/T2 -4000*maxinEtOH/T2 -4000*maxinEtOH/T2/3 ...
                      0 omega0-J J-0.01 -pi -pi -pi -40*maxinEtOH -2000*maxinEtOH];
                ub = [4000*maxinEtOH/T2 4000*maxinEtOH/T2 4000*maxinEtOH/T2/3 ...
                      T2*100 omega0+J J+0.01 pi pi pi 40*maxinEtOH 1000*maxinEtOH];
                
                % Least-squares model fitting
                ThreeLorentzModelInit = lsqcurvefit(@ThreeLorentzModel, ThreeLorentzModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)), lb, ub, lsqopts);
                [ThreeLorentzModelParam, resid] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)), @ThreeLorentzModel, ThreeLorentzModelInit, nlinopts);
                
                EtOHheight = max(ThreeLorentzModelParam(1:3)) * ThreeLorentzModelParam(4);
                MRS_struct.out.(vox{kk}).(target{trg}).FitError(ii) = 100*std(resid)/EtOHheight;
                MRS_struct.out.(vox{kk}).(target{trg}).Area(ii) = sum(ThreeLorentzModel(ThreeLorentzModelParam,freq(freqbounds))) * abs(freq(1) - freq(2));
                MRS_struct.out.(vox{kk}).(target{trg}).FWHM(ii) = 1./(pi*ThreeLorentzModelParam(4))*1e3;
                MRS_struct.out.(vox{kk}).(target{trg}).ModelParam(ii,:) = ThreeLorentzModelParam;
                MRS_struct.out.(vox{kk}).(target{trg}).Resid(ii,:) = resid;
                
                % Calculate SNR of EtOH signal (MM: 170502)
                noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                MRS_struct.out.(vox{kk}).(target{trg}).SNR(ii) = abs(EtOHheight)/noiseSigma_DIFF;
                
            else
                
                error('Fitting MRS_struct.p.target not recognised');
                
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   1a. Start up the output figure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ishandle(102)
                clf(102); % MM (170629)
            end
            h = figure(102);
            % MM (170629): Open figure in center of screen
            scr_sz = get(0, 'ScreenSize');
            fig_w = 1000;
            fig_h = 707;
            set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
            set(h,'Color',[1 1 1]);
            figTitle = 'GannetFit Output';
            set(gcf,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
            
            % Spectra plot
            ha = subplot(2,2,1);
            metabmin = min(real(DIFF(ii,plotbounds)));
            metabmax = max(real(DIFF(ii,plotbounds)));
            resmax = max(resid);
            resid = resid + metabmin - resmax;
            if strcmp(target{trg},'GABA')
                plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                    freq(freqbounds), ThreeLorentzModel(ThreeLorentzModelParam,freq(freqbounds)), 'r', ...
                    freq(freqbounds), resid, 'k');
                set(gca,'XLim',[2.6 3.4]);
            elseif strcmp(target{trg},'GSH')
                plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b' ,...
                    freq(freqbounds), SixLorentzModel(SixLorentzModelParam,freq(freqbounds)), 'r', ...
                    freq(freqbounds),resid, 'k');
                set(gca,'XLim',[2.5 3.4]);
            elseif strcmp(target{trg},'Lac')
                plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                    freq(freqbounds), FourGaussModel(FourGaussModelParam,freq(freqbounds)), 'r', ...
                    freq(freqbounds), FourGaussModel(MMGaussModelParam,freq(freqbounds)), 'r' , ...
                    freq(freqbounds), resid, 'k');
                set(gca,'XLim',[0 2.1]);
            elseif strcmp(target{trg},'Glx')
                plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                    freq(freqbounds), TwoLorentzModel2(TwoLorentzModelParam,freq(freqbounds)), 'r', ...
                    freq(freqbounds), resid, 'k');
                set(gca,'XLim',[3.35 4.15]);
            elseif strcmp(target{trg},'EtOH')
                plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                    freq(freqbounds), ThreeLorentzModel(ThreeLorentzModelParam,freq(freqbounds)), 'r', ...
                    freq(freqbounds), resid, 'k');
                set(gca,'XLim',[0.7 1.6]);
            end
                        
            title('Edited Spectrum and Model Fit');
            set(gca,'XDir','reverse');
            
            % From here on is cosmetic - adding labels etc.
            switch target{trg}
                case 'GABA'
                    h1 = text(3.15,metabmax/3,MRS_struct.p.target);
                    set(h1, 'horizontalAlignment', 'center');
                    labelbounds = freq <= 2.8 & freq >= 2.6; % MM (170705)
                    tailtop = max(real(DIFF(ii,labelbounds)));
                    tailbottom = min(real(DIFF(ii,labelbounds)));
                    h2 = text(2.775, min(resid), 'residual');
                    set(h2, 'horizontalAlignment', 'left');
                    text(2.775, tailtop+metabmax/20, 'data', 'Color', [0 0 1]);
                    text(2.775, tailbottom-metabmax/20, 'model', 'Color', [1 0 0]);
                    
                case 'GSH'
                    h1 = text(3.05,maxinGSH/2,target{trg});
                    set(h1, 'horizontalAlignment', 'center');
                    labelbounds = freq <= 2.8 & freq >= 2.5; % MM (170705)
                    tailtop = max(real(DIFF(ii,labelbounds)));
                    tailbottom = min(real(DIFF(ii,labelbounds)));
                    h2 = text(2.7, min(resid),'residual');
                    set(h2, 'horizontalAlignment', 'left');
                    text(2.7, tailtop+metabmax/20, 'data', 'Color', [0 0 1]);
                    text(2.7, tailbottom-20*metabmax/20, 'model', 'Color', [1 0 0]);
                    
                case 'Glx'
                    h1 = text(MRS_struct.out.(vox{kk}).(target{trg}).ModelParam(ii,4), metabmax/5, MRS_struct.p.target);
                    set(h1, 'horizontalAlignment', 'center');
                    labelbounds = freq <= 3.6 & freq >= 3.4; % MM (170705)
                    tailtop = max(real(DIFF(ii,labelbounds)));
                    tailbottom = min(real(DIFF(ii,labelbounds)));
                    h2 = text(3.5, min(resid),'residual');
                    set(h2, 'horizontalAlignment', 'left');
                    text(3.5, tailtop+metabmax/20, 'data', 'Color', [0 0 1]);
                    text(3.5, tailbottom-metabmax/20, 'model', 'Color', [1 0 0]);
                    
                case 'EtOH'
                    h1 = text(1.35,metabmax/3,MRS_struct.p.target);
                    set(h1, 'horizontalAlignment', 'center');
                    labelbounds = freq <= 1.0 & freq >= 0.7; % MM (170705)
                    tailtop = max(real(DIFF(ii,labelbounds)));
                    tailbottom = min(real(DIFF(ii,labelbounds)));
                    h2 = text(0.9, min(resid), 'residual');
                    set(h2, 'horizontalAlignment', 'left');
                    text(0.9, tailtop+metabmax/20, 'data', 'Color', [0 0 1]);
                    text(0.9, tailbottom-metabmax/20, 'model', 'Color', [1 0 0]);
            end
            xlabel('ppm');
            set(gca,'YTick',[]);
            set(gca,'Box','off');
            set(gca,'YColor','white');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   2.  Water Fit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp(MRS_struct.p.Reference_compound,'H2O')
                                
                % Estimate height and baseline from data
                [maxinWater, watermaxindex] = max(real(WaterData(ii,:)),[],2);
                waterbase = mean(real(WaterData(ii,freq <= 4 & freq >= 3.8)));
                                
                LGPModelInit = [maxinWater 20 freq(watermaxindex) 0 waterbase -50 0];
                lb = [0.01*maxinWater 1 4.6 0 0 -50 -pi];
                ub = [40*maxinWater 100 5.0 0.000001 1 0 pi];
                
                freqbounds = freq <= 5.1 & freq >= 4.5;
                plotbounds = freq <= 5.2 & freq >= 4.4;
                
                % Least-squares model fitting
                LGPModelInit = lsqcurvefit(@LorentzGaussModelP, LGPModelInit, freq(freqbounds), real(WaterData(ii,freqbounds)), lb, ub, lsqopts);
                [LGPModelParam, residw] = nlinfit(freq(freqbounds), real(WaterData(ii,freqbounds)), @LorentzGaussModelP, LGPModelInit, nlinopts);
                
                WaterArea = sum(real(LorentzGaussModel(LGPModelParam(1:end-1),freq(freqbounds))) - BaselineModel(LGPModelParam(3:5),freq(freqbounds)),2);
                MRS_struct.out.(vox{kk}).water.Area(ii) = WaterArea * abs(freq(1)-freq(2));
                waterheight = LGPModelParam(1);
                MRS_struct.out.(vox{kk}).water.FitError(ii) = 100*std(residw)/waterheight;
                % MM (170202)
                LG = real(LorentzGaussModel(LGPModelParam(1:end-1),freq(freqbounds))) - BaselineModel(LGPModelParam(3:5),freq(freqbounds));
                LG = LG./max(LG);
                ind = find(LG >= 0.5);
                f = freq(freqbounds);
                w = abs(f(ind(1)) - f(ind(end)));
                MRS_struct.out.(vox{kk}).water.FWHM(ii) = w * MRS_struct.p.LarmorFreq(ii);
                MRS_struct.out.(vox{kk}).water.ModelParam(ii,:) = LGPModelParam;
                MRS_struct.out.(vox{kk}).water.Resid(ii,:) = residw; % MM (160913)
                
                % Calculate SNR of water signal (MM: 170502)
                noiseSigma_Water = CalcNoise(freq, WaterData(ii,:));
                MRS_struct.out.(vox{kk}).water.SNR(ii) = abs(waterheight)/noiseSigma_Water;
                
                % Root sum square fit error and concentration in institutional units -- MGSaleh & MM
                switch target{trg}
                    case 'GABA'
                        MRS_struct.out.(vox{kk}).GABA.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).GABA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                        MRS_struct = CalcConc(MRS_struct, vox{kk}, 'GABA', ii);
                        
                    case 'Glx'
                        MRS_struct.out.(vox{kk}).Glx.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).Glx.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                        MRS_struct = CalcConc(MRS_struct, vox{kk}, 'Glx', ii);
                                                
                    case 'GSH'
                        MRS_struct.out.(vox{kk}).GSH.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).GSH.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                        MRS_struct = CalcConc(MRS_struct, vox{kk}, (target{trg}), ii);
                        
                    case 'Lac'
                        MRS_struct.out.(vox{kk}).Lac.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).Lac.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                        MRS_struct = CalcConc(MRS_struct, vox{kk}, (target{trg}), ii);
                        
                    case 'EtOH'
                        MRS_struct.out.(vox{kk}).EtOH.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).EtOH.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                        MRS_struct = CalcConc(MRS_struct, vox{kk}, 'EtOH', ii);
                end
                
                % Generate scaled spectra (for plotting) CJE Jan2011, MM (170705)
                MRS_struct.spec.(vox{kk}).(target{trg}).off_scaled(ii,:) = ...
                    MRS_struct.spec.(vox{kk}).(target{trg}).off(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                MRS_struct.spec.(vox{kk}).(target{trg}).on_scaled(ii,:) = ...
                    MRS_struct.spec.(vox{kk}).(target{trg}).on(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                MRS_struct.spec.(vox{kk}).(target{trg}).diff_scaled(ii,:) = ...
                    MRS_struct.spec.(vox{kk}).(target{trg}).diff(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                
                % MM (170703): Reorder structure fields 
                MRS_struct.out.(vox{kk}).water = orderfields(MRS_struct.out.(vox{kk}).water, {'Area', 'FWHM', 'SNR', 'ModelParam', 'Resid', 'FitError'});
                
                hb = subplot(2,2,3);
                watmin = min(real(WaterData(ii,:)));
                watmax = max(real(WaterData(ii,:)));
                resmax = max(residw);
                residw = residw + watmin - resmax;
                plot(freq(plotbounds), real(WaterData(ii,plotbounds)), 'b', ...
                    freq(freqbounds), real(LorentzGaussModelP(LGPModelParam,freq(freqbounds))), 'r', ...
                    freq(freqbounds), residw, 'k');
                set(gca,'XDir','reverse');
                set(gca,'YTick',[]);
                set(gca,'Box','off');
                set(gca,'YColor','white');
                xlim([4.4 5.2]);
                % Add on some labels
                hwat = text(4.85,watmax/2,'Water');
                set(hwat,'HorizontalAlignment','right');
                % Get the right vertical offset for the residual label
                labelfreq = freq(freqbounds); % MM (170705)
                rlabelbounds = labelfreq <= 4.7 & labelfreq >= 4.4;
                axis_bottom = axis;
                hwatres = text(4.6, max(min(residw(rlabelbounds))-0.05*watmax, axis_bottom(3)), 'residual');
                set(hwatres, 'horizontalAlignment', 'left');
                
            end
            
            % MM (170703): Reorder structure fields
            if ~MRS_struct.p.HERMES % MM (170703): work on this for GSH data
                if strcmp(MRS_struct.p.Reference_compound,'H2O')
                    MRS_struct.out.(vox{kk}).(target{trg}) = orderfields(MRS_struct.out.(vox{kk}).(target{trg}), ...
                        {'Area', 'FWHM', 'SNR', 'ModelParam', 'Resid', 'FitError', 'FitError_W', 'ConcIU'});
                else
                    MRS_struct.out.(vox{kk}).(target{trg}) = orderfields(MRS_struct.out.(vox{kk}).(target{trg}), ...
                        {'Area', 'FWHM', 'SNR', 'ModelParam', 'Resid', 'FitError'});
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   5. Build GannetFit Output
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp(MRS_struct.p.Reference_compound,'H2O')
                set(gca,'Box','off');
                set(gca,'YColor','white');
                xlabel('ppm');
                title('Reference Signal');
            end
            
            % And running the plot
            if any(strcmp('mask',fieldnames(MRS_struct))) == 1
                hc = subplot(2,2,2);
                get(hc,'pos'); % get position of axes
                set(hc,'pos',[0.52 0.52 0.42 0.42]) % move the axes slightly
                size_max = size(MRS_struct.mask.img{ii},1);
                imagesc(MRS_struct.mask.img{ii}(:,size_max+(1:size_max)));
                colormap('gray');
                caxis([0 1])
                axis equal;
                axis tight;
                axis off;
                subplot(2,2,4,'replace');
            else
                subplot(2,2,2);
                axis off;
            end
            
            % MM (170703): Cleaner text alignment
            text_pos = 0.9; % A variable to determine y-position of text on printout on figure -- Added by MGSaleh
                        
            % MM (180112)
            if strcmp(MRS_struct.p.vendor,'Siemens_rda')
                [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii*2-1});
            else
                [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii});
            end
            
            text(0, text_pos, 'Filename', 'FontName', 'Helvetica', 'FontSize', 10);
            text(0.375, text_pos, [': ' tmp tmp2], 'FontName', 'Helvetica', 'FontSize', 10, 'Interpreter', 'none');
            
            % Some changes to accomodate multiplexed fitting output
            switch target{trg}
                case 'GABA'
                    tmp1 = 'GABA+ Area';
                    tmp2 = sprintf(': %.3g', MRS_struct.out.(vox{kk}).GABA.Area(ii));
                    
                case 'Glx'
                    tmp1 = 'Glx Area';
                    tmp2 = sprintf(': %.3g', MRS_struct.out.(vox{kk}).Glx.Area(ii));
                      
                case 'GSH'
                    tmp1 = 'GSH Area';
                    tmp2 = sprintf(': %.3g', MRS_struct.out.(vox{kk}).(target{trg}).Area(ii));
                    
                case 'Lac'
                    tmp1 = 'Lac Area';
                    tmp2 = sprintf(': %.3g', MRS_struct.out.(vox{kk}).(target{trg}).Area(ii));
                    
                case 'EtOH'
                    tmp1 = 'EtOH Area';
                    tmp2 = sprintf(': %.3g', MRS_struct.out.(vox{kk}).EtOH.Area(ii));
            end
            text(0, text_pos-0.1, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
            text(0.375, text_pos-0.1, tmp2, 'FontName', 'Helvetica', 'FontSize', 10);
            
            if strcmp(MRS_struct.p.Reference_compound,'H2O')
                tmp = sprintf(': %.3g', MRS_struct.out.(vox{kk}).water.Area(ii));
                text(0, text_pos-0.2, 'Water Area', 'FontName', 'Helvetica', 'FontSize', 10);
                text(0.375, text_pos-0.2, tmp, 'FontName', 'Helvetica', 'FontSize', 10);
                
                tmp = sprintf(': %.1f Hz', MRS_struct.out.(vox{kk}).water.FWHM(ii));
                text(0, text_pos-0.3, 'Water FWHM', 'FontName', 'Helvetica', 'FontSize', 10);
                text(0.375, text_pos-0.3, tmp, 'FontName', 'Helvetica', 'FontSize', 10);
                
                tmp1 = sprintf(': %.1f%%', MRS_struct.out.(vox{kk}).(target{trg}).FitError_W(ii));
                tmp2 = [target{trg} ' Conc.'];
                tmp3 = sprintf(': %.3f mM', MRS_struct.out.(vox{kk}).(target{trg}).ConcIU(ii));
                
                text(0, text_pos-0.4, 'Fit Error', 'FontName', 'Helvetica', 'FontSize', 10);
                text(0.375, text_pos-0.4, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
                text(0, text_pos-0.5, tmp2, 'FontName', 'Helvetica', 'FontSize', 10);
                text(0.375, text_pos-0.5, tmp3, 'FontName', 'Helvetica', 'FontSize', 10);
            else
                tmp = sprintf(': %.1f Hz', MRS_struct.out.(vox{kk}).(target{trg}).FWHM(ii));
                text(0, text_pos-0.3, [target{trg} ' FWHM'], 'FontName', 'Helvetica', 'FontSize', 10);
                text(0.375, text_pos-0.3, tmp, 'FontName', 'Helvetica', 'FontSize', 10);
                
                tmp1 = sprintf(': %.1f%%', MRS_struct.out.(vox{kk}).(target{trg}).FitError(ii));
                text(0, text_pos-0.4, 'Fit Error', 'FontName', 'Helvetica', 'FontSize', 10);
                text(0.375, text_pos-0.4, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
            end
            
            tmp = [': ' MRS_struct.version.fit];
            text(0, text_pos-0.6, 'FitVer', 'FontName', 'Helvetica', 'FontSize', 10);
            text(0.375, text_pos-0.6, tmp, 'FontName', 'Helvetica', 'FontSize', 10);
            
            % Add Gannet logo
            if any(strcmp('mask',fieldnames(MRS_struct)))
                subplot(2,2,4);
            else
                subplot(2,2,4,'replace');
            end
            axis off;
            script_path = which('GannetFit');
            Gannet_logo = [script_path(1:(end-12)) '/Gannet3_logo.png'];
            A2 = imread(Gannet_logo,'png','BackgroundColor',[1 1 1]);
            axes('Position',[0.80, 0.05, 0.15, 0.15]);
            image(A2);
            axis off;
            axis square;            

            % For Philips .data
            if strcmpi(MRS_struct.p.vendor,'Philips_data')
                fullpath = MRS_struct.metabfile{ii};
                fullpath = regexprep(fullpath, '.data', '_data');
                fullpath = regexprep(fullpath, '\', '_');
                fullpath = regexprep(fullpath, '/', '_');
            end
            
            % MM (180112)
            if strcmp(MRS_struct.p.vendor,'Siemens_rda')
                [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii*2-1});
            else
                [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii});
            end            
            
            if sum(strcmp(listfonts,'Helvetica')) > 0
                if strcmp(MRS_struct.p.Reference_compound,'H2O')
                    set([ha,hb],'FontName','Helvetica');
                else
                    set(ha,'FontName','Helvetica');
                end
            end
            
            % Save PDF output
            set(gcf,'PaperUnits','inches');
            set(gcf,'PaperSize',[11 8.5]);
            set(gcf,'PaperPosition',[0 0 11 8.5]);
            
            if ~exist('GannetFit_output','dir')
                mkdir GannetFit_output;
            end
            
            if strcmpi(MRS_struct.p.vendor,'Philips_data')
                pdfname = fullfile('GannetFit_output', [fullpath '_' target{trg} '_' vox{kk} '_fit.pdf']); % MM (180112)
            else
                pdfname = fullfile('GannetFit_output', [metabfile_nopath '_' target{trg} '_' vox{kk} '_fit.pdf']); % MM (180112)
            end            
            saveas(h, pdfname);
            
        end
        
    end
    
    % 140116: ADH reorder structure
    if isfield(MRS_struct, 'mask')
        if isfield(MRS_struct, 'waterfile')
            structorder = {'version', 'ii', ...
                'metabfile', 'waterfile', 'p', 'fids', 'spec', 'out', 'mask'};
        else
            structorder = {'version', 'ii', ...
                'metabfile', 'p', 'fids', 'spec', 'out', 'mask'};
        end
    else
        if isfield(MRS_struct, 'waterfile')
            structorder = {'version', 'ii', ...
                'metabfile', 'waterfile', 'p', 'fids', 'spec', 'out'};
        else
            structorder = {'version','ii', ...
                'metabfile', 'p', 'fids', 'spec', 'out'};
        end
    end
    MRS_struct = orderfields(MRS_struct, structorder);
    
    if MRS_struct.p.mat % save MRS_struct as mat file
        mat_name = ['MRS_struct_' vox{kk} '.mat'];
        save(mat_name,'MRS_struct');
    end
    
    if MRS_struct.p.csv % export MRS_struct fields into csv file
        ExportToCSV(MRS_struct, target, kk, 'fit');
    end
    
    fprintf('\n\n');
    
end


%%%%%%%%%%%%%%%%%%%% THREE LORENTZ MODEL (MM: 180723) %%%%%%%%%%%%%%%%%%%%
function F = ThreeLorentzModel(x,freq)
% ThreeLorentzModel with phase parameters
% Based on Marshall & Roe, 1978 (Analytical Chem)

Ha    = x(1);  % amplitude of outer peak
Hb    = x(2);  % amplitude of outer peak
H0    = x(3);  % amplitude of middle peak
T2    = x(4);  % T2 relaxation time constant
f0    = x(5);  % frequency of middle peak (in ppm)
J     = x(6);  % J-coupling constant (in ppm)
phi_a = x(7);  % phase of outer peak (in rad)
phi_b = x(8);  % phase of outer peak
phi_0 = x(9);  % phase of middle peak
M     = x(10); % baseline slope
C     = x(11); % baseline offset

Aa = cos(phi_a) .* ((Ha .* T2) ./ (1 + (f0 + J - freq).^2 .* T2.^2)) ...
    - sin(phi_a) .* ((Ha .* (f0 + J - freq) .* T2.^2) ./ (1 + (f0 + J - freq).^2 .* T2.^2));

Ab = cos(phi_b) .* ((Hb .* T2) ./ (1 + (f0 - J - freq).^2 .* T2.^2)) ...
    - sin(phi_b) .* ((Hb .* (f0 - J - freq) .* T2.^2) ./ (1 + (f0 - J - freq).^2 .* T2.^2));

A0 = cos(phi_0) .* ((H0 .* T2) ./ (1 + (f0 - freq).^2 .* T2.^2)) ...
    - sin(phi_0) .* ((H0 .* (f0 - freq) .* T2.^2) ./ (1 + (f0 - freq).^2 .* T2.^2));

F = Aa + A0 + Ab + M .* (f0 - freq) + C;


%%%%%%%%%%%%%%%%%%%% TWO LORENTZ MODEL (MM: 180723) %%%%%%%%%%%%%%%%%%%%
function F = TwoLorentzModel2(x,freq)
% TwoLorentzModel with phase parameters
% Based on Marshall & Roe, 1978 (Analytical Chem)

Ha    = x(1); % amplitude of peak 1
Hb    = x(2); % amplitude of peak 2
T2    = x(3); % T2 relaxation time constant
f0    = x(4); % frequency (in ppm)
J     = x(5); % J-coupling constant (in ppm)
phi_a = x(6); % phase of peak 1 (in rad)
phi_b = x(7); % phase of peak 2
M     = x(8); % baseline slope
C     = x(9); % baseline offset

Aa = cos(phi_a) .* ((Ha .* T2) ./ (1 + (f0 + J - freq).^2 .* T2.^2)) ...
    - sin(phi_a) .* ((Ha .* (f0 + J - freq) .* T2.^2) ./ (1 + (f0 + J - freq).^2 .* T2.^2));

Ab = cos(phi_b) .* ((Hb .* T2) ./ (1 + (f0 - J - freq).^2 .* T2.^2)) ...
    - sin(phi_b) .* ((Hb .* (f0 - J - freq) .* T2.^2) ./ ( 1 + (f0 - J - freq).^2 .* T2.^2));

F = Aa + Ab + M .* (f0 - freq) + C;


%%%%%%%%%%%%%%%%%%%% SIX LORENTZ MODEL (MM: 180912) %%%%%%%%%%%%%%%%%%%%
function F = SixLorentzModel(x,freq)
% TwoLorentzModel with phase parameters
% Based on Marshall & Roe, 1978 (Analytical Chem)

Ha    = x(1); % amplitude of peak 1
Hb    = x(2); % amplitude of peak 2
Hc    = x(3); % amplitude of peak 3
Hd    = x(4); % amplitude of peak 4
He    = x(5); % amplitude of peak 5
Hf    = x(6); % amplitude of peak 6
T2    = x(7); % T2 relaxation time constant
f0    = x(8); % frequency (in ppm)
J     = x(9); % J-coupling constant (in ppm)
phi_a = x(10); % phase of peak 1 (in rad)
phi_b = x(11); % phase of peak 2
phi_c = x(12); % phase of peak 1 (in rad)
phi_d = x(13); % phase of peak 2
phi_e = x(14); % phase of peak 1 (in rad)
phi_f = x(15); % phase of peak 2
M     = x(16); % baseline slope
C     = x(17); % baseline offset

Aa = cos(phi_a) .* ((Ha .* T2) ./ (1 + (f0 + 3*J - freq).^2 .* T2.^2)) ...
    - sin(phi_a) .* ((Ha .* (f0 + 3*J - freq) .* T2.^2) ./ (1 + (f0 + 3*J - freq).^2 .* T2.^2));

Ab = cos(phi_b) .* ((Hb .* T2) ./ (1 + (f0 + 2*J - freq).^2 .* T2.^2)) ...
    - sin(phi_b) .* ((Hb .* (f0 + 2*J - freq) .* T2.^2) ./ ( 1 + (f0 + 2*J - freq).^2 .* T2.^2));


Ac = cos(phi_c) .* ((Hc .* T2) ./ (1 + (f0 + J - freq).^2 .* T2.^2)) ...
    - sin(phi_c) .* ((Hc .* (f0 + J - freq) .* T2.^2) ./ (1 + (f0 + J - freq).^2 .* T2.^2));

Ad = cos(phi_d) .* ((Hd .* T2) ./ (1 + (f0 - J - freq).^2 .* T2.^2)) ...
    - sin(phi_d) .* ((Hd .* (f0 - J - freq) .* T2.^2) ./ ( 1 + (f0 - J - freq).^2 .* T2.^2));


Ae = cos(phi_e) .* ((He .* T2) ./ (1 + (f0 - 2*J - freq).^2 .* T2.^2)) ...
    - sin(phi_e) .* ((He .* (f0 - 2*J - freq) .* T2.^2) ./ (1 + (f0 - 2*J - freq).^2 .* T2.^2));

Af = cos(phi_f) .* ((Hf .* T2) ./ (1 + (f0 - 3*J - freq).^2 .* T2.^2)) ...
    - sin(phi_f) .* ((Hf .* (f0 - 3*J - freq) .* T2.^2) ./ ( 1 + (f0 - 3*J - freq).^2 .* T2.^2));

F = Aa + Ab + Ac + Ad + Ae + Af + M .* (f0 - freq) + C;


%%%%%%%%%%%%%%%%  LORENTZGAUSSMODEL %%%%%%%%%%%%%%%%%%%%
function F = LorentzGaussModel(x,freq)
% Function for LorentzGaussModel Model

% CJE 24Nov10 - removed phase term from fit - this is now dealt with
% by the phasing of the water ref scans in MRSLoadPfiles
%Lorentzian Model multiplied by a Gaussian.
% x(1) = Amplitude of (scaled) Lorentzian
% x(2) = 1 / hwhm of Lorentzian (hwhm = half width at half max)
% x(3) = centre freq of Lorentzian
% x(4) = linear baseline slope
% x(5) = constant baseline amplitude
% x(6) =  -1 / 2 * sigma^2  of gaussian

% Lorentzian  = (1/pi) * (hwhm) / (deltaf^2 + hwhm^2)
% Peak height of Lorentzian = 4 / (pi*hwhm)
% F is a normalised Lorentzian - height independent of hwhm
%   = Lorentzian / Peak

F = (x(1)*ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)) ... % Lorentzian
    .* (exp(x(6)*(freq-x(3)).*(freq-x(3)))) ... % Gaussian
    + x(4)*(freq-x(3)) ... % linear baseline
    + x(5); % constant baseline


%%%%%%%%%%%%%%%% LORENTZGAUSSMODEL WITH PHASE %%%%%%%%%%%%%%%%%%%%
function F = LorentzGaussModelP(x,freq)
% Function for LorentzGaussModel Model with Phase

% Lorentzian Model multiplied by a Gaussian
% x(1) = Amplitude of (scaled) Lorentzian
% x(2) = 1 / hwhm of Lorentzian (hwhm = half width at half max)
% x(3) = centre freq of Lorentzian
% x(4) = linear baseline slope
% x(5) = constant baseline amplitude
% x(6) =  -1 / 2 * sigma^2  of gaussian
% x(7) = phase (in rad)

% Lorentzian  = (1/pi) * (hwhm) / (deltaf^2 + hwhm^2)
% Peak height of Lorentzian = 4 / (pi*hwhm)
% F is a normalised Lorentzian - height independent of hwhm
%   = Lorentzian / Peak

F = ((cos(x(7))*x(1)*ones(size(freq)) + sin(x(7))*x(1)*x(2)*(freq-x(3)))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)) ... % Lorentzian
    .* (exp(x(6)*(freq-x(3)).*(freq-x(3)))) ... % Gaussian
    + x(4)*(freq-x(3)) ... % linear baseline
    + x(5); % constant baseline


%%%%%%%%%%%%%%% BASELINE %%%%%%%%%%%%%%%%%%%%%%%
function F = BaselineModel(x,freq)
% Function for Baseline Model

F = x(2)*(freq-x(1))+x(3);


%%%%%%%%%%%%%%%%%%% CALC CONC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MRS_struct = CalcConc(MRS_struct, vox, metab, ii)
% Function for quantifying concentration in absolute units
% Convert metabolits and water areas to institutional units
% (concentration in mmol/L)

TR = MRS_struct.p.TR(ii)/1000;
TE = MRS_struct.p.TE(ii)/1000;
if isfield(MRS_struct.p,'TR_water')
    TR_water = MRS_struct.p.TR_water(ii)/1000;
else
    TR_water = TR;
end
if isfield(MRS_struct.p,'TE_water')
    TE_water = MRS_struct.p.TE_water(ii)/1000;
else
    TE_water = TE;
end
PureWaterConc = 55000; % mmol/L
T1_Water = 2928/1e3; % from unpublished JHU data (K. Chan)
T2_Water = 0.503; % NB: T2 of water in CSF in vivo; Piechnik et al. 2009 (MRM)
N_H_Water = 2;

switch metab
    case 'GABA'
        EditingEfficiency = 0.5; % for TE = 68 ms
        T1_Metab = 1.84; % Harris et al. 2017 (MRI)
        T2_Metab = 248/1e3; % Harris et al. 2017 (MRI)
        N_H_Metab = 2;
        
    case 'Glx'
        EditingEfficiency = 0.4; % determined by FID-A simulations (for TE = 68 ms)
        T1_Metab = 1.18; % Choi et al. 2006 (MRM)
        T2_Metab = 640/1e3; % Choi et al. 2006 (MRM)
        N_H_Metab = 1;
        
    case 'GSH'
        EditingEfficiency = 0.74;  % At 3T based on Quantification of Glutathione in the Human Brain by MR Spectroscopy at 3 Tesla:
        % Comparison of PRESS and MEGA-PRESS
        % Faezeh Sanaei Nezhad etal. DOI 10.1002/mrm.26532, 2016 -- MGSaleh
        T1_Metab = 0.40; % At 3T based on Doubly selective multiple quantum chemical shift imaging and
        % T1 relaxation time measurement of glutathione (GSH) in the human brain in vivo
        % In-Young Choi et al. NMR Biomed. 2013; 26: 28?34 -- MGSaleh
        T2_Metab = 0.12; % At 3T based on the ISMRM abstract
        % T2 relaxation times of 18 brain metabolites determined in 83 healthy volunteers in vivo
        % Milan Scheidegger et al. Proc. Intl. Soc. Mag. Reson. Med. 22 (2014)-- MGSaleh
        N_H_Metab = 2;
        
    case 'Lac'
        EditingEfficiency = 0.94; % determined by FID-A simulations (for TE = 140 ms)
        T1_Metab = 1.50; % Wijnen et al. 2015 (NMR Biomed)
        T2_Metab = 0.24; % Madan et al. 2015 (MRM) (NB: this was estimated in brain tumors)
        N_H_Metab = 3;
        
    case 'EtOH'
        EditingEfficiency = 0.5; % determined by FID-A simulations (for TE = 140 ms)
        T1_Metab = 1.50; % Wijnen et al. 2015 (NMR Biomed)
        T2_Metab = 0.24; % Madan et al. 2015 (MRM) (NB: this was estimated in brain tumors)
        N_H_Metab = 2;
end

T1_Factor = (1-exp(-TR_water./T1_Water)) ./ (1-exp(-TR./T1_Metab));
T2_Factor = exp(-TE_water./T2_Water) ./ exp(-TE./T2_Metab);

if strcmpi(MRS_struct.p.vendor,'Siemens_rda')
    % Factor of 2 is appropriate for averaged Siemens data (read in separately as ON and OFF)
    MRS_struct.out.(vox).(metab).ConcIU(ii) = (MRS_struct.out.(vox).(metab).Area(ii) ./ MRS_struct.out.(vox).water.Area(ii))  ...
        .* PureWaterConc .* T1_Factor .* T2_Factor .* (N_H_Water./N_H_Metab) ...
        ./ 2 ./ EditingEfficiency;
else
    MRS_struct.out.(vox).(metab).ConcIU(ii) = (MRS_struct.out.(vox).(metab).Area(ii) ./ MRS_struct.out.(vox).water.Area(ii))  ...
        .* PureWaterConc .* T1_Factor .* T2_Factor .* (N_H_Water./N_H_Metab) ...
        ./ EditingEfficiency;
end



