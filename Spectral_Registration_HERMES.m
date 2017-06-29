function [AllFramesFTrealign, MRS_struct] = Spectral_Registration_HERMES(MRS_struct)
% Spectral registration is a time-domain frequency-and-phase correction
% routine as per Near et al. (2015). Incorporates a multiplexed,
% probabilistic approach for aligning HERMES data (MM: 170609)

showPlots = 'n';

% Looping parameters
if MRS_struct.p.HERMES % run registration four times - once for each HERMES experiment
    SpecRegLoop = 3;
    SubspecToAlign = repmat([3 2 1 0], [1 size(MRS_struct.fids.data,2)/4]);
else % run registration once or twice for MEGA-PRESS acquisitions
    %SpecRegLoop = 1;
    %SubspecToAlign = repmat([1 0], [1 size(MRS_struct.fids.data,2)/2]);
    SpecRegLoop = 0;
    SubspecToAlign = zeros(1, size(MRS_struct.fids.data,2));
end

% Pre-allocate memory
ii = MRS_struct.ii;
MRS_struct.fids.data_align = zeros(size(MRS_struct.fids.data));
MRS_struct.out.SpecReg.ParsFit = zeros(size(MRS_struct.fids.data,2),2);
zMSE = zeros(1,size(MRS_struct.fids.data,2));
CorrParsML = zeros(size(MRS_struct.fids.data,2),2);
count = 0;
parsGuess = [0 0];

% Inputs
DataToAlign = MRS_struct.fids.data;
time = (0:1:(MRS_struct.p.npoints-1)).'/MRS_struct.p.sw;
input.dwelltime = 1/MRS_struct.p.sw;

% Probability density function and parameter bounds
Cauchy = @(x,s,l) s./(pi.*(s.^2+(x-l).^2));
lb = [0 -Inf];
ub = [Inf Inf];

% Optimization options
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts, 'MaxIter', 1e4);
mleopts  = statset('MaxIter', 1e4, 'MaxFunEvals', 1e4);

% Set dimensions of figures of histograms
if strcmpi(showPlots,'y')
    d.w = 0.6;
    d.h = 0.45;
    d.l = (1-d.w)/2;
    d.b = (1-d.h)/2;
end

while SpecRegLoop > -1
    
    % Use first N points of time-domain data, where N is the last point where SNR > 3
    noise = std(real(DataToAlign(ceil(0.75*size(DataToAlign(:,SubspecToAlign == SpecRegLoop),1)):end,SubspecToAlign == SpecRegLoop)),[],1);
    noise = mean(noise);
    signal = mean(abs(DataToAlign(:,SubspecToAlign == SpecRegLoop)),2);
    SNR = signal./noise;
    N = find(SNR > 3);
    tMax = N(end);
    
    % 'Flatten' complex data for use in nlinfit
    clear flatdata;
    flatdata(:,1,:) = real(DataToAlign(1:tMax,SubspecToAlign == SpecRegLoop));
    flatdata(:,2,:) = imag(DataToAlign(1:tMax,SubspecToAlign == SpecRegLoop));
    
    % Select reference transient
    flattarget = median(flatdata,3); % median across transients
    %flattarget = squeeze(flatdata(:,:,1)); % Mth transient
    target = flattarget(:);
    
    % Pre-allocate memory
    if ~count
        parsFit = zeros(size(flatdata,3), 2);
        MSE = zeros(1, size(flatdata,3));
    end
    
    % Determine frequency and phase offsets by spectral registration
    reverseStr = '';
    for corrloop = 1:size(flatdata,3)
        msg = sprintf('\nSpectral registration - Fitting transient: %d', corrloop);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
        transient = squeeze(flatdata(:,:,corrloop));
        input.data = transient(:);
        [parsFit(corrloop,:), ~, ~, ~, MSE(corrloop)] = nlinfit(input, target, @FreqPhaseShiftNest, parsGuess, nlinopts);
        parsGuess = parsFit(corrloop,:);
    end
    
    count = count + 1;
    
    % Probability distribution of frequency offsets (estimated by maximum likelihood)
    MRS_struct.out.MLalign.f.x(count,:,ii) = parsFit(:,1);
    start = [std(MRS_struct.out.MLalign.f.x(count,:,ii))/2, median(MRS_struct.out.MLalign.f.x(count,:,ii))];
    [MRS_struct.out.MLalign.f.p(count,:,ii), MRS_struct.out.MLalign.f.p_ci(:,:,count,ii)] = ...
        mle(MRS_struct.out.MLalign.f.x(count,:,ii), 'pdf', Cauchy, 'start', start, 'lower', lb, 'upper', ub, 'options', mleopts);
    MRS_struct.out.MLalign.f.fx(count,:,ii) = ...
        linspace(1.5*min(MRS_struct.out.MLalign.f.x(count,:,ii)), 1.5*max(MRS_struct.out.MLalign.f.x(count,:,ii)), 1e3);
    MRS_struct.out.MLalign.f.pdf(count,:,ii) = Cauchy(MRS_struct.out.MLalign.f.fx(count,:,ii), ...
        MRS_struct.out.MLalign.f.p(count,1,ii), MRS_struct.out.MLalign.f.p(count,2,ii));
    
    % Probability distribution of phase offsets (estimated by maximum likelihood)
    MRS_struct.out.MLalign.ph.x(count,:,ii) = parsFit(:,2);
    start = [std(MRS_struct.out.MLalign.ph.x(count,:,ii))/2, median(MRS_struct.out.MLalign.ph.x(count,:,ii))];
    [MRS_struct.out.MLalign.ph.p(count,:,ii), MRS_struct.out.MLalign.ph.p_ci(:,:,count,ii)] = ...
        mle(MRS_struct.out.MLalign.ph.x(count,:,ii), 'pdf', Cauchy, 'start', start, 'lower', lb, 'upper', ub, 'options', mleopts);
    MRS_struct.out.MLalign.ph.fx(count,:,ii) = ...
        linspace(1.5*min(MRS_struct.out.MLalign.ph.x(count,:,ii)), 1.5*max(MRS_struct.out.MLalign.ph.x(count,:,ii)), 1e3);
    MRS_struct.out.MLalign.ph.pdf(count,:,ii) = Cauchy(MRS_struct.out.MLalign.ph.fx(count,:,ii), ...
        MRS_struct.out.MLalign.ph.p(count,1,ii), MRS_struct.out.MLalign.ph.p(count,2,ii));
    
    if strcmpi(showPlots,'y')
        % Histogram of frequency offsets
        H = figure(333);
        set(H, 'Color', 'w', 'Units', 'Normalized', 'OuterPosition', [d.l d.b d.w d.h]);
        subplot(1,2,1);
        bins = linspace(min(MRS_struct.out.MLalign.f.x(count,:,ii)), max(MRS_struct.out.MLalign.f.x(count,:,ii)), 15);
        binWidth = abs(bins(1)-bins(2));
        h = bar(bins, histc(MRS_struct.out.MLalign.f.x(count,:,ii),bins)/(length(MRS_struct.out.MLalign.f.x(count,:,ii))*binWidth), 'histc');
        h.FaceColor = [0.8 0.8 0.8];
        hold on;
        plot(MRS_struct.out.MLalign.f.fx(count,:,ii), MRS_struct.out.MLalign.f.pdf(count,:,ii), 'Color', [1 0 0], 'LineWidth', 1.2);
        hold off;
        xlabel('\Deltaf (Hz)', 'FontSize', 15);
        ylabel('P(\theta)', 'FontSize', 15);
        set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
        
        % Histogram of phase offsets
        subplot(1,2,2);
        bins = linspace(min(MRS_struct.out.MLalign.ph.x(count,:,ii)), max(MRS_struct.out.MLalign.ph.x(count,:,ii)), 15);
        binWidth = abs(bins(1)-bins(2));
        h = bar(bins, histc(MRS_struct.out.MLalign.ph.x(count,:,ii),bins)/(length(MRS_struct.out.MLalign.ph.x(count,:,ii))*binWidth), 'histc');
        h.FaceColor = [0.8 0.8 0.8];
        hold on
        plot(MRS_struct.out.MLalign.ph.fx(count,:,ii), MRS_struct.out.MLalign.ph.pdf(count,:,ii), 'Color', [1 0 0], 'LineWidth', 1.2)
        hold off
        xlabel('\Delta\phi (deg)', 'FontSize', 15);
        set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
        
        drawnow;
        pause(1);
    end
    
    corrloop_d = find(SubspecToAlign == SpecRegLoop);
    MRS_struct.out.SpecReg.ParsFit(corrloop_d,:,ii) = parsFit;
    CorrParsML(corrloop_d,1) = parsFit(:,1) - MRS_struct.out.MLalign.f.p(count,2,ii)';
    CorrParsML(corrloop_d,2) = parsFit(:,2) - MRS_struct.out.MLalign.ph.p(count,2,ii)';
    zMSE(corrloop_d) = zscore(MSE); % standardized MSEs
    
    % Apply frequency and phase corrections
    for corrloop = 1:size(flatdata,3)
        % Default correction
        %DataToAlign(:,corrloop_d(corrloop)) = DataToAlign(:,corrloop_d(corrloop)) .* ...
        %    exp(1i*parsFit(corrloop,1)*2*pi*time) * exp(1i*pi/180*parsFit(corrloop,2));
        
        % Freq/phase correction + Cauchy pdf location parameter shift
        DataToAlign(:,corrloop_d(corrloop)) = DataToAlign(:,corrloop_d(corrloop)) .* ...
            exp(1i*(parsFit(corrloop,1) - MRS_struct.out.MLalign.f.p(count,2,ii))*2*pi*time) * ...
            exp(1i*pi/180*(parsFit(corrloop,2) - MRS_struct.out.MLalign.ph.p(count,2,ii)));
    end
    
    if SpecRegLoop == 0
        
%         MRS_struct.out.MLalign.f_aligned.x(ii,:) = CorrParsML(:,1);
%         start = [std(MRS_struct.out.MLalign.f_aligned.x(ii,:))/2, median(MRS_struct.out.MLalign.f_aligned.x(ii,:))];
%         MRS_struct.out.MLalign.f_aligned.p(ii,:) = ...
%             mle(MRS_struct.out.MLalign.f_aligned.x(ii,:), 'pdf', Cauchy, 'start', start, 'lower', lb.c, 'upper', ub.c, 'options', mleopts);
%         MRS_struct.out.MLalign.f_aligned.fx(ii,:) = ...
%             linspace(1.1*min(MRS_struct.out.MLalign.f_aligned.x(ii,:)), 1.1*max(MRS_struct.out.MLalign.f_aligned.x(ii,:)), 1e3);
%         MRS_struct.out.MLalign.f_aligned.pdf(ii,:) = ...
%             Cauchy(MRS_struct.out.MLalign.f_aligned.fx(ii,:), MRS_struct.out.MLalign.f_aligned.p(ii,1), MRS_struct.out.MLalign.f_aligned.p(ii,2));
%         
%         MRS_struct.out.MLalign.ph_aligned.x(ii,:) = CorrParsML(:,2);
%         start = [std(MRS_struct.out.MLalign.ph_aligned.x(ii,:))/2, median(MRS_struct.out.MLalign.ph_aligned.x(ii,:))];
%         MRS_struct.out.MLalign.ph_aligned.p(ii,:) = ...
%             mle(MRS_struct.out.MLalign.ph_aligned.x(ii,:), 'pdf', Cauchy, 'start', start, 'lower', lb.c, 'upper', ub.c, 'options', mleopts);
%         MRS_struct.out.MLalign.ph_aligned.fx(ii,:) = ...
%             linspace(1.1*min(MRS_struct.out.MLalign.ph_aligned.x(ii,:)), 1.1*max(MRS_struct.out.MLalign.ph_aligned.x(ii,:)), 1e3);
%         MRS_struct.out.MLalign.ph_aligned.pdf(ii,:) = ...
%             Cauchy(MRS_struct.out.MLalign.ph_aligned.fx(ii,:), MRS_struct.out.MLalign.ph_aligned.p(ii,1), MRS_struct.out.MLalign.ph_aligned.p(ii,2));
%         
%         clf(H);
%         subplot(1,2,1);
%         bins = linspace(min(MRS_struct.out.MLalign.f_aligned.x(ii,:)), max(MRS_struct.out.MLalign.f_aligned.x(ii,:)), 20);
%         binWidth = abs(bins(1)-bins(2));
%         h = bar(bins, histc(MRS_struct.out.MLalign.f_aligned.x(ii,:),bins)/(length(MRS_struct.out.MLalign.f_aligned.x(ii,:))*binWidth), 'histc');
%         h.FaceColor = [0.8 0.8 0.8];
%         hold on;
%         plot(MRS_struct.out.MLalign.f_aligned.fx(ii,:), MRS_struct.out.MLalign.f_aligned.pdf(ii,:), 'Color', [1 0 0], 'LineWidth', 1.2);
%         hold off;
%         xlabel('\Deltaf (Hz)', 'FontSize', 15);
%         ylabel('P(x)', 'FontSize', 15);
%         set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
%         
%         subplot(1,2,2);
%         bins = linspace(min(MRS_struct.out.MLalign.ph_aligned.x(ii,:)), max(MRS_struct.out.MLalign.ph_aligned.x(ii,:)), 20);
%         binWidth = abs(bins(1)-bins(2));
%         h = bar(bins, histc(MRS_struct.out.MLalign.ph_aligned.x(ii,:),bins)/(length(MRS_struct.out.MLalign.ph_aligned.x(ii,:))*binWidth), 'histc');
%         h.FaceColor = [0.8 0.8 0.8];
%         hold on;
%         plot(MRS_struct.out.MLalign.ph_aligned.fx(ii,:), MRS_struct.out.MLalign.ph_aligned.pdf(ii,:), 'Color', [1 0 0], 'LineWidth', 1.2);
%         hold off;
%         xlabel('\Delta\phi (deg)', 'FontSize', 15);
%         ylabel('P(x)', 'FontSize', 15);
%         set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
%         
%         drawnow;
%         %pause(1);        
      
        % Line-broadening, zero-filling and FFT
        FullData = DataToAlign .* repmat((exp(-time*MRS_struct.p.LB*pi)), [1 size(MRS_struct.fids.data,2)]);
        AllFramesFTrealign = fftshift(fft(FullData,MRS_struct.p.ZeroFillTo,1),1);
        
        % In frequency domain, shift Cr signals to 3.02 and get frequency 'right' as opposed to 'consistent'
        freqrange = MRS_struct.spec.freq >= 2.925 & MRS_struct.spec.freq <= 3.125;
        %freqrange = MRS_struct.spec.freq >= 3.125 & MRS_struct.spec.freq <= 3.32; % Cho
        
        [~,FrameMaxPos] = max(real(AllFramesFTrealign(freqrange,:)),[],1);
        freq = MRS_struct.spec.freq(freqrange);
        CrFreqShift = freq(FrameMaxPos);
        CrFreqShift = CrFreqShift - 3.02; % 3.2
        CrFreqShift_pts = round(CrFreqShift / abs(MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2)));
        
        % Apply circular frequency shifts
        for corrloop = 1:size(AllFramesFTrealign,2)
            AllFramesFTrealign(:,corrloop) = circshift(AllFramesFTrealign(:,corrloop), CrFreqShift_pts(corrloop));
        end        
        
        % Reject transients that are greater than +/-3 st. devs. of MSEs
        rejectFrames = zMSE > 3 | zMSE < -3;
        if MRS_struct.p.HERMES
            rejectFrames = reshape(rejectFrames, [4 size(MRS_struct.fids.data,2)/4])';
        else
            rejectFrames = reshape(rejectFrames, [2 size(MRS_struct.fids.data,2)/2])';
        end
        ind = any(rejectFrames,2);
        rejectFrames(ind,:) = 1;
        MRS_struct.out.reject(:,ii) = reshape(rejectFrames', [1 size(MRS_struct.fids.data,2)]);
        
    end
    
    SpecRegLoop = SpecRegLoop - 1;
    
end

if exist('H','var')
    close(H);
end
fprintf('\n');

end


% function output = FreqPhaseShiftNest(pars,input)
% % Heavily copied from jn_freqPhaseShiftNest of Jamie Near, McGill.
% 
% f = pars(1); % frequency shift [Hz]
% p = pars(2); % phase shift [deg]
% 
% dwelltime = input.dwelltime;
% t = 0:dwelltime:(length(input.data)/2-1)*dwelltime;
% input.data = reshape(input.data, [length(input.data)/2, 2]);
% fid = input.data(:,1) + 1i*input.data(:,2);
% 
% y = fid.*exp(1i*pi*(t'*f*2+p/180));
% output = [real(y), imag(y)];
% output = output(:);
% 
% end

