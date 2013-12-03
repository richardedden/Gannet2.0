    function [FitParams, rejectframe, residCr]  = FitPeaksByFrames2(freq, FrameData, initx)
    
    %options = optimset('lsqcurvefit');
    %options = optimset(options,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',1e5 , ...
    %		   'MaxFunEvals', 1e12);
    % initx = [area hwhm f0 phase baseline0 baseline1]
    
    nlinopts = statset('nlinfit');
    nlinopts = statset(nlinopts, 'MaxIter', 1e4, 'Display','Off');
    nframes = size(FrameData,2);
    FrameData=real(FrameData);
    for jj = 1:nframes
        % [fit_param, resnorm, resid, exitflag ]  = ...
        %     lsqcurvefit(@(xdummy,ydummy) LorentzModel(xdummy, ydummy), initx, ...
        % 		  freq', real(FrameData(:,jj)));
        [fit_param, residCr] = nlinfit(freq', (FrameData(:,jj)), ...
            @(xdummy, ydummy) LorentzModel(xdummy, ydummy), ...
            initx, nlinopts);
        FitParams(jj,:) = fit_param;
        FitParams2(jj,:) = fit_param;
        fit_plot = LorentzModel(fit_param, freq);
        %fit_plot = LorentzModel_nolinear(fit_param, freq);
        
          %figure(3); subplot(1,2,1);plot(freq', real(FrameData(:,jj)), 'g', freq', fit_plot,'b');
          %pause(2)
        %set(gca,'XDir','reverse');
        %  input('next')
    end
    
    for kk=1:size(FitParams,1)
        if FitParams(kk,1)<0
            FitParams(kk,4)= FitParams(kk,4)+pi;
        end
    end
    
    % Need to deal with phase wrap:
    % Convert to complex number then recalculate phase within 2*pi range
    phase_wrapped = FitParams(:,4);
    cmplx = cos(phase_wrapped) + 1i * sin(phase_wrapped);
    phase_unwrapped = angle(cmplx);
    

    % then fix to be within -pi..pi
    offsetpos =  2*pi*lt(phase_unwrapped, -pi);
    offsetneg = -2*pi*gt(phase_unwrapped,  pi);
    phase_unwrapped = phase_unwrapped + offsetpos + offsetneg;
    FitParams(:,4) = phase_unwrapped;
    
    % Fix area and linewidth to be positive

    FitParams(:,1) = abs(FitParams(:,1));
    FitParams(:,2) = abs(FitParams(:,2));
    
    %for jj = 1:nframes
    %    fitplot1=LorentzModel(FitParams2(jj,:), freq);;
    %    figure(3); subplot(1,2,1);plot(freq', real(FrameData(:,jj)), 'g', freq', fitplot1,'b');
        %fitplot2=ComplexLorentzModel(FitParams(jj,:), freq);;
        %figure(3); subplot(1,2,2);plot(freq', real(FrameData(:,jj)), 'g', freq', fitplot2,'b');pause(2);
    %end
    % Reject any point where the fit params - area, fwhm, phase
    %  or freq are > 3stdev away from the mean
    % set reject criteria for all fit parameters
    MeanFitParams = mean(FitParams, 1);
    UpperLim = repmat(MeanFitParams + 3*std(FitParams,1), [nframes 1]);
    LowerLim = repmat(MeanFitParams - 3*std(FitParams,1), [nframes 1]);
    %but don't reject on linear, const baseline fit vals
    UpperLim(:,5:6) = Inf;
    LowerLim(:,5:6) = -Inf;
    rejectframe = gt(FitParams, UpperLim);
    rejectframe = rejectframe + lt(FitParams, LowerLim);
    rejectframe = max(rejectframe,[],2);
    
    end
