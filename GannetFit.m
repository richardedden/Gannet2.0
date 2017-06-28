function [MRS_struct] = GannetFit(MRS_struct, varargin)
%
% MRS_struct = structure with data loaded from MRSLoadPfiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gannet 3.0 version of Gannet Fit - analysis tool for GABA-edited MRS.
% Need some new sections like
%   1. GABA, Glx and GSH Fit -- MGSaleh January 2017
%   2. Water Fit
%   3. Cr Fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Changes made to the code for Gannet3.0:
% GABA and Glx outputs are saved in seperate mat structures -- MGSaleh 29 June 2016
% Concentration estimates from GABAGlx fitting -- MGSaleh 13 July 2016
% GSH output is saved in seperate mat structures -- MGSaleh January 2017


%% Input parameters and some definitions:
%Determine whether multiple regions and fitting targets or not and then use it to
%create a loop -- MGSaleh 2016

if MRS_struct.p.PRIAM    % Deciding how many regions are there -- MGSaleh 2016 
    
    
    reg={MRS_struct.p.Reg};
    
else
    reg={MRS_struct.p.Reg{1}};
    
end


if (MRS_struct.p.HERMES) & nargin < 2
    target  ={MRS_struct.p.target, MRS_struct.p.target2};
    varargin={MRS_struct.p.target, MRS_struct.p.target2};
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
            case 'GABAGlx'
                MRS_struct.p.target = 'GABAGlx';
            case 'GSH'
                MRS_struct.p.target = 'GSH';
        end
    end
    
    target = {MRS_struct.p.target};
end





FIT_LSQCURV = 0;
FIT_NLINFIT = 1;
fit_method = FIT_NLINFIT; %FIT_NLINFIT;
waterfit_method = FIT_NLINFIT;

% if strcmp(MRS_struct.p.target, 'GABAGlx')
%     
%     Datadiff = MRS_struct.spec.GABAGlx.diff; % Added by MGSaleh 2016
%     Dataoff  = MRS_struct.spec.GABAGlx.off; % Added by MGSaleh 2016
%     Dataon   = MRS_struct.spec.GABAGlx.on; % Added by MGSaleh 2016
%     
% else
%     
%     Datadiff = MRS_struct.spec.GSH.diff; % Added by MGSaleh 2016
%     Dataoff  = MRS_struct.spec.GSH.off; % Added by MGSaleh 2016
%     Dataon   = MRS_struct.spec.GSH.on; % Added by MGSaleh 2016
%     
%     
% end





freq=MRS_struct.spec.freq;
if strcmp(MRS_struct.p.Reference_compound,'H2O')
    WaterData=MRS_struct.spec.water;
end
MRS_struct.versionfit = '161029'; %format - yy/mm/dd  -- MGSaleh 2016
disp(['GABA Fit Version is ' MRS_struct.versionfit ]);
fitwater=1;
% eval(['MRS_struct.spec.', target{1}, '.diff(ii,:)'])
% numscans=size(Datadiff);
% numscans=numscans(1);
numscans=size(MRS_struct.spec.(MRS_struct.p.Reg{1}).(target{1}).diff); % Added by MGSaleh 2016
numscans=numscans(1);
%110624
epsdirname = [ './MRSfit_' datestr(clock,'yymmdd') ];


%% Metabolite fitting 
for kk = 1:length(reg)
    
for trg=1:length(target)
        
        %Defining variables -- MGSaleh 2016
        Datadiff = MRS_struct.spec.(reg{kk}).(target{trg}).diff; % Added by MGSaleh 2016
        Dataoff  = MRS_struct.spec.(reg{kk}).(target{trg}).off;  % Added by MGSaleh 2016 % e.g: MRS_struct.spec.GABAGlx.off; % Added by MGSaleh 2016
        Dataon   = MRS_struct.spec.(reg{kk}).(target{trg}).on;   % Added by MGSaleh 2016 % e.g: MRS_struct.spec.GABAGlx.on; % Added by MGSaleh 2016
        
        numscans=size(Datadiff);
        numscans=numscans(1);
       
        disp(sprintf('\n\n Working on: %s',target{trg}));
        
        for ii=1:numscans
        MRS_struct.gabafile{ii};
       
    

    if strcmp(target{trg},'GABA')
    
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % 1.  GABA Fit 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        % ...from GaussModel;
        % x(1) = gaussian amplitude
        % x(2) = 1/(2*sigma^2)
        % x(3) = centre freq of peak
        % x(4) = amplitude of linear baseline
        % x(5) = constant amplitude offset

        %Hard code it to fit from 2.79 ppm to 3.55 ppm
        z=abs(MRS_struct.spec.freq-3.55);
        lowerbound=find(min(z)==z);
        z=abs(MRS_struct.spec.freq-2.79);%2.75
        upperbound=find(min(z)==z);
        freqbounds=lowerbound:upperbound;
        plotbounds=(lowerbound-150):(upperbound+150);
        maxinGABA=abs(max(real(Datadiff(MRS_struct.ii,freqbounds)))-min(real(Datadiff(MRS_struct.ii,freqbounds))));
        
        % smarter estimation of baseline params, Krish's idea (taken from Johns
        % code; NAP 121211
        grad_points = (real(Datadiff(ii,upperbound)) - real(Datadiff(ii,lowerbound))) ./ ...
            (upperbound - lowerbound); %in points
        LinearInit = grad_points ./ (MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2)); %in ppm
        constInit = (real(Datadiff(ii,upperbound)) + real(Datadiff(ii,lowerbound))) ./2;
        xval = [ 1:(upperbound-lowerbound+1) ];
        linearmodel = grad_points .* xval + Datadiff(ii,lowerbound);
        %End copy code
        resnorm=zeros([numscans size(freqbounds,2)]);
        GaussModelInit = [maxinGABA -90 3.026 -LinearInit constInit]; %default in 131016
        lb = [0 -200 2.87 -40*maxinGABA -2000*maxinGABA]; %NP; our bounds are 0.03 less due to creatine shift
        ub = [4000*maxinGABA -40 3.12 40*maxinGABA 1000*maxinGABA];
        %plot(freq(freqbounds),real(Datadiff(ii,freqbounds)),freq(freqbounds),GaussModel_area(GaussModelInit,freq(freqbounds)))
        options = optimset('lsqcurvefit');
        options = optimset(options,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',1e5);
        nlinopts = statset('nlinfit');
        nlinopts = statset(nlinopts, 'MaxIter', 1e5);
         
        %Fitting to a Gaussian model happens here
         [GaussModelParam(ii,:),resnorm,residg] = lsqcurvefit(@(xdummy,ydummy) GaussModel_area(xdummy,ydummy), ...
            GaussModelInit, freq(freqbounds),real(Datadiff(ii,freqbounds)), ...
            lb,ub,options);
            residg = -residg;
        if(fit_method == FIT_NLINFIT)
            GaussModelInit = GaussModelParam(ii,:);
            % 1111013 restart the optimisation, to ensure convergence
            for fit_iter = 1:100
                [GaussModelParam(ii,:), residg, J, COVB, MSE] = nlinfit(freq(freqbounds), real(Datadiff(ii,freqbounds)), ... % J, COBV, MSE edited in
                    @(xdummy,ydummy) GaussModel_area(xdummy,ydummy), ...
                    GaussModelInit, ...
                    nlinopts);
                MRS_struct.out.GABA.fitparams_iter(fit_iter,:,ii) = GaussModelParam(ii,:); %specifically for GABA - MGSaleh
                GaussModelInit = GaussModelParam(ii,:);
                ci = nlparci(GaussModelParam(ii,:), residg,'covar',COVB); %copied over
            end
        end
        GABAheight = GaussModelParam(ii,1);
        % FitSTD reports the standard deviation of the residuals / gaba HEIGHT
        MRS_struct.out.(reg{kk}).(target{trg}).FitError(ii)  =  100*std(residg)/GABAheight;
        % This sets GabaArea as the area under the curve.
        MRS_struct.out.(reg{kk}).(target{trg}).Area(ii)=GaussModelParam(ii,1)./sqrt(-GaussModelParam(ii,2))*sqrt(pi);
        sigma = ( 1 / (2 * (abs(GaussModelParam(ii,2)))) ).^(1/2);
        MRS_struct.out.(reg{kk}).(target{trg}).FWHM(ii) =  abs( (2* MRS_struct.p.LarmorFreq) * sigma);
        MRS_struct.out.(reg{kk}).(target{trg}).ModelFit(ii,:)=GaussModelParam(ii,:);
        MRS_struct.out.(reg{kk}).(target{trg}).resid(ii,:) = residg;
        MRS_struct.out.(reg{kk}).(target{trg}).snr(ii) = GABAheight / std(residg);  % Added by MGSaleh

        
        
    elseif strcmp (target{trg},'GSH') %strcmp(MRS_struct.p.target,'GSH')
        %GSH fit
        
        %Hard code it to fit from 2 ppm to 4 ppm
        z=abs(MRS_struct.spec.freq-3.3);
        lowerbound=find(min(z)==z);
        z=abs(MRS_struct.spec.freq-2.35);
        upperbound=find(min(z)==z);
        freqbounds=lowerbound:upperbound;
        plotbounds=(lowerbound-150):(upperbound+150);
        offset = real(mean(Datadiff(ii, freqbounds(1:10)),2) + mean(Datadiff(ii, freqbounds((end-9):end)),2))/2;
        slope = real(mean(Datadiff(ii, freqbounds(1:10)),2) - mean(Datadiff(ii, freqbounds((end-9):end)),2))/real(MRS_struct.spec.freq(freqbounds(1)) - MRS_struct.spec.freq(freqbounds(end)));

        peak_amp = 0.03; %Presumably this won't work for some data... for now it seems to work.

        initx = [peak_amp*0.7 -300 2.95 peak_amp*0.8 -500 2.73 -peak_amp*2.5 -1000 2.61 -peak_amp*2.5 -1000 2.55 peak_amp*0.5 -600 2.42  0 0 0.02 ];
        lb = [0 -5000   2.9 0 -5000  2.68 -0.3 -5000    2.57 -0.3 -5000   2.48 0 -5000    2.3 -0.1 -0.01 -0.01 ];
        ub = [0.1 -50   3.0 0.1 -50   2.8   0 -50    2.68      0 -50   2.57 0.1 0     2.43 0.1 0.01 0.03 ]; 
        figure(87)
        plot(MRS_struct.spec.freq(freqbounds),FiveGaussModel(initx,MRS_struct.spec.freq(freqbounds)),MRS_struct.spec.freq(freqbounds),real(Datadiff(ii,freqbounds)));%MRS_struct.spec.freq(freqbounds),real(Datadiff(ii,freqbounds)
        options = optimset('lsqcurvefit');
        options = optimset(options,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',1e5);
        nlinopts = statset('nlinfit');
        nlinopts = statset(nlinopts, 'MaxIter', 1e5);
        
        [FiveGaussModelParam(ii,:),resnorm,residg] = lsqcurvefit(@(xdummy,ydummy) FiveGaussModel(xdummy,ydummy), ...
        initx, MRS_struct.spec.freq(freqbounds),real(Datadiff(ii,freqbounds)), ...
        lb,ub,options);
        residg = -residg;
        
%         if(fit_method == FIT_NLINFIT)
%             initx = FiveGaussModelParam(ii,:);
%             % 1111013 restart the optimisation, to ensure convergence
%             for fit_iter = 1:10
%                 [FiveGaussModelParam(ii,:),resnorm,residg] = lsqcurvefit(@(xdummy,ydummy) FiveGaussModel(xdummy,ydummy), ...
%                     initx, MRS_struct.spec.freq(freqbounds),real(Datadiff(ii,freqbounds)), ...
%                     lb,ub,options);
% %                 MRS_struct.out.GABA.fitparams_iter(fit_iter,:,ii) = FiveGaussModelParams(ii,:); %specifically for GABA - MGSaleh
%                 initx = FiveGaussModelParam(ii,:);
% %                 ci = nlparci(GaussModeFiveGaussModelParamlParam(ii,:), residg,'covar',COVB); %copied over
%             end
%         end
        
        
        MRS_struct.out.(reg{kk}).GSH.FiveGaussModel(ii,:)=FiveGaussModelParam(ii,:);
        GSHGaussModelParam(ii,:)=FiveGaussModelParam(ii,:);
        GSHGaussModelParam(ii,4:3:13)=0;
        NAAGaussModelParam(ii,:)=FiveGaussModelParam(ii,:);
        NAAGaussModelParam(1)=0;

        BaselineModelParam=GSHGaussModelParam;
        BaselineModelParam(ii,1)=0;

%         MRS_struct.out.(reg{kk}).(target{trg}).Area(ii)=real(sum(FiveGaussModel(GSHGaussModelParam(ii,:), MRS_struct.spec.freq(freqbounds))-FiveGaussModel(BaselineModelParam(ii,:), MRS_struct.spec.freq(freqbounds))))*(MRS_struct.spec.freq(1)-MRS_struct.spec.freq(2));
%         %Not sure how to handle fit error. For now, do whole range
%         GSHheight = GSHGaussModelParam(ii,1);
%         MRS_struct.out.(reg{kk}).(target{trg}).FitError(ii)=  100*std(residg)/GSHheight;
%         
        % GSH fitting output -- MGSaleh
        MRS_struct.out.(reg{kk}).(target{trg}).Area(ii) = real(sum(FiveGaussModel(GSHGaussModelParam(ii,:), MRS_struct.spec.freq(freqbounds))-FiveGaussModel(BaselineModelParam(ii,:), MRS_struct.spec.freq(freqbounds))))*(MRS_struct.spec.freq(1)-MRS_struct.spec.freq(2));
        GSHheight = GSHGaussModelParam(ii,1);
        MRS_struct.out.(reg{kk}).(target{trg}).Height(ii) = GSHheight;
        
        % Range to determine residuals for GSH 
        z=abs(MRS_struct.spec.freq-3.3);  % For GSH - Added by MGSaleh 2017
        midbound1=find(min(z)==z);        % For GSH - Added by MGSaleh 2017
        z=abs(MRS_struct.spec.freq-2.77);
        midbound2=find(min(z)==z);
        
        MRS_struct.out.(reg{kk}).(target{trg}).FitError(ii)  =  100*std(residg(:,[1:(midbound2-(midbound1-1))]))/GSHheight; % The residuals is from the whole range specified above for the the fitting --MGSaleh
        sigma = ( 1 / (2 * (abs(GSHGaussModelParam(ii,2)))) ).^(1/2);
        MRS_struct.out.(reg{kk}).(target{trg}).FWHM(ii) =  abs( (2* MRS_struct.p.LarmorFreq) * sigma);
        MRS_struct.out.(reg{kk}).(target{trg}).ModelFit(ii,:)=GSHGaussModelParam(ii,:); % GaussModelParam(ii,7:9) are related to GABA  -- MGSaleh;
        MRS_struct.out.(reg{kk}).(target{trg}).resid(ii,:) = residg(:,[1:(midbound2-(midbound1-1))]);
        MRS_struct.out.(reg{kk}).(target{trg}).snr(ii) = GSHheight / std(residg(:,[1:(midbound2-(midbound1-1))]));  % Added by MGSaleh

    
    
    
        
        
    elseif strcmp (target{trg},'Lac') % strcmp(MRS_struct.p.target,'Lac')    
        %This is a work in progress - currenly mainly code copied form GSH
        %Hard code it to fit from 0.8 ppm to 1.8 ppm
        z=abs(MRS_struct.spec.freq-1.8);
        lowerbound=find(min(z)==z);
        z=abs(MRS_struct.spec.freq-0.5);
        upperbound=find(min(z)==z);
        freqbounds=lowerbound:upperbound;
        plotbounds=(lowerbound-150):(upperbound+150);
        offset = (mean(Datadiff(ii, freqbounds(1:10)),2) + mean(Datadiff(ii, freqbounds((end-9):end)),2))/2;
        slope = (mean(Datadiff(ii, freqbounds(1:10)),2) - mean(Datadiff(ii, freqbounds((end-9):end)),2))/(MRS_struct.spec.freq(freqbounds(1)) - MRS_struct.spec.freq(freqbounds(end)));
        peak_amp = 0.03; %Presumably this won't work for some data... for now it seems to work.

        initx = [peak_amp*0.16 -100 1.18 peak_amp*0.3 -1000 1.325   offset slope 0];
        lb = [0 -300 0.9 0 -5000 1.0  -1 -1 -1];
        ub = [1 0 1.4 1 0 1.6  1 1 1]; 
        
        options = optimset('lsqcurvefit');
        options = optimset(options,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',1e5);
        nlinopts = statset('nlinfit');
        nlinopts = statset(nlinopts, 'MaxIter', 1e5);
        
        %Plot the starting fit model
        %figure(99)
        %subplot(1,2,1)
        %plot(MRS_struct.spec.freq(freqbounds),real(MRS_struct.spec.diff(ii,freqbounds)),MRS_struct.spec.freq(freqbounds),FourGaussModel(initx,MRS_struct.spec.freq(freqbounds)));

        
        
        [FourGaussModelParam(ii,:),resnorm,residg] = lsqcurvefit(@(xdummy,ydummy) FourGaussModel(xdummy,ydummy), ...
        initx, MRS_struct.spec.freq(freqbounds),real(Datadiff(ii,freqbounds)), ...
        lb,ub,options);
    
        %subplot(1,2,2)
        %plot(MRS_struct.spec.freq(freqbounds),real(MRS_struct.spec.diff(ii,freqbounds)),MRS_struct.spec.freq(freqbounds),FourGaussModel(FourGaussModelParam,MRS_struct.spec.freq(freqbounds)),MRS_struct.spec.freq(freqbounds),FourGaussModel([FourGaussModelParam(1:3) 0 FourGaussModelParam(5:end)],MRS_struct.spec.freq(freqbounds)));
        %error('Fitting Init model plot')
    
    
        residg = -residg;
        LacGaussModelParam(ii,:)=FourGaussModelParam(ii,:);
        LacGaussModelParam(ii,1)=0;
        MMGaussModelParam(ii,:)=FourGaussModelParam(ii,:);
        MMGaussModelParam(4)=0;

        BaselineModelParam=MMGaussModelParam;
        BaselineModelParam(ii,1)=0;

        MRS_struct.out.GABAArea(ii)=real(sum(FourGaussModel(LacGaussModelParam(ii,:), MRS_struct.spec.freq(freqbounds))-FourGaussModel(BaselineModelParam(ii,:), MRS_struct.spec.freq(freqbounds))))*(MRS_struct.spec.freq(1)-MRS_struct.spec.freq(2));
        %Not sure how to handle fit error. For now, do whole range
        GABAheight = LacGaussModelParam(ii,4);
        MRS_struct.out.GABAFitError(ii)=  100*std(residg)/GABAheight;
        
       
    
    elseif strcmp (MRS_struct.p.target,'Glx')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Glx Fit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Hard code it to fit from 3.45 ppm to 4.10 ppm
        % MM: Larger fit range helps to avoid fitting of Cr-CH3 artefact(?)
        % at 3.91 ppm
        z=abs(MRS_struct.spec.freq-4.10); %4.10
        lowerbound=find(min(z)==z);
        z=abs(MRS_struct.spec.freq-3.45); %3.45
        upperbound=find(min(z)==z);
        freqbounds=lowerbound:upperbound;
        plotbounds=(lowerbound-150):(upperbound+150);
        maxinGABA=max(real(Datadiff(MRS_struct.ii,freqbounds)));
        % smarter estimation of baseline params, Krish's idea (taken from Johns
        % code; NAP 121211
        grad_points = (real(Datadiff(ii,upperbound)) - real(Datadiff(ii,lowerbound))) ./ ...
            (upperbound - lowerbound); %in points
        LinearInit = grad_points ./ (MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2)); %in ppm
        constInit = (real(Datadiff(ii,upperbound)) + real(Datadiff(ii,lowerbound))) ./2;
        xval = [ 1:(upperbound-lowerbound+1) ];
        linearmodel = grad_points .* xval + Datadiff(ii,lowerbound);
        %End copy code
        resnorm=zeros([numscans size(freqbounds,2)]);
        
        % To fit a Double Gaussian

            % MM: Allowing peaks to vary individually seems to work better
            % than keeping the distance fixed (i.e., including J in the
            % function)
            GaussModelInit = [maxinGABA -90 3.72 maxinGABA -90 3.77 -LinearInit constInit];
            lb = [0 -200 3.72-0.01 0 -200 3.77-0.01 -40*maxinGABA -2000*maxinGABA];
            ub = [4000*maxinGABA -40 3.72+0.01 4000*maxinGABA -40 3.77+0.01 40*maxinGABA 1000*maxinGABA];

        
        options = optimset('lsqcurvefit');
        options = optimset(options,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',1e5);
        nlinopts = statset('nlinfit');
        nlinopts = statset(nlinopts, 'MaxIter', 1e5);
        
        
            [GaussModelParam(ii,:),resnorm,residg] = lsqcurvefit(@(xdummy,ydummy) DoubleGaussModel_area(xdummy,ydummy), ...
                GaussModelInit, freq(freqbounds),real(Datadiff(ii,freqbounds)), ...
                lb,ub,options);
            residg = -residg;
            if(fit_method == FIT_NLINFIT)
                GaussModelInit = GaussModelParam(ii,:);
                % 111013 restart the optimisation, to ensure convergence
                for fit_iter = 1:100
                    [GaussModelParam(ii,:), residg, J, COVB, MSE] = nlinfit(freq(freqbounds), real(Datadiff(ii,freqbounds)), ... % J, COBV, MSE edited in
                        @(xdummy,ydummy) DoubleGaussModel_area(xdummy,ydummy), ...
                        GaussModelInit, ...
                        nlinopts);
                    MRS_struct.out.Glx.fitparams_iter(fit_iter,:,ii) = GaussModelParam(ii,:);
                    GaussModelInit = GaussModelParam(ii,:);
                    ci = nlparci(GaussModelParam(ii,:), residg,'covar',COVB); %copied over
                end
            end

            
            Glxheight = max(GaussModelParam(ii,[1,4]));
            % GABAFitError reports the standard deviation of the residuals / GABAheight
            MRS_struct.out.(reg{kk}).(target{trg}).FitError(ii) = 100*std(residg)/Glxheight;
            % This sets GABAArea as the area under the curve
            MRS_struct.out.(reg{kk}).(target{trg}).Area(ii) = (GaussModelParam(ii,1)./sqrt(-GaussModelParam(ii,2))*sqrt(pi)) + ...
                (GaussModelParam(ii,4)./sqrt(-GaussModelParam(ii,5))*sqrt(pi));
            sigma = ((1/(2*(abs(GaussModelParam(ii,2))))).^(1/2)) + ((1/(2*(abs(GaussModelParam(ii,5))))).^(1/2));
            MRS_struct.out.(reg{kk}).(target{trg}).FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq)*sigma);
            MRS_struct.out.(reg{kk}).(target{trg}).ModelFit(ii,:) = GaussModelParam(ii,:);
            MRS_struct.out.(reg{kk}).(target{trg}).resid(ii,:) = residg;
            MRS_struct.out.(reg{kk}).(target{trg}).snr(ii) = Glxheight / std(residg);
        
        
    elseif strcmp (target{trg},'GABAGlx') %strcmp (MRS_struct.p.target,'GABAGlx')  % Need to quantify GABA and Glx separately  --MGSaleh
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  TEsting a GABA and Glx Fit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Hard code it to fit from 2.79 ppm to 4.10 ppm
        % Full range of GABA and Glx fitting above
        z=abs(MRS_struct.spec.freq-4.05); %4.10
        lowerbound=find(min(z)==z);
        z=abs(MRS_struct.spec.freq-2.79); %2.79
        upperbound=find(min(z)==z);
        freqbounds=lowerbound:upperbound;
        plotbounds=(lowerbound-150):(upperbound+150);
        
        % Range to determine residuals for GABA and Glx separately
        z=abs(MRS_struct.spec.freq-3.55); % For GABA - Added by MGSaleh 2016
        midbound1=find(min(z)==z);        % For GABA - Added by MGSaleh 2016
        
        z=abs(MRS_struct.spec.freq-3.45); % For Glx  - Added by MGSaleh 2016
        midbound2=find(min(z)==z);        % For Glx  - Added by MGSaleh 2016
        
        
        maxinGABA=max(real(Datadiff(MRS_struct.ii,freqbounds)));
        % smarter estimation of baseline params, Krish's idea (taken from Johns
        % code; NAP 121211
        grad_points = (real(Datadiff(ii,upperbound)) - real(Datadiff(ii,lowerbound))) ./ ...
            (upperbound - lowerbound); %in points
        LinearInit = grad_points ./ (MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2)); %in ppm
        constInit = (real(Datadiff(ii,upperbound)) + real(Datadiff(ii,lowerbound))) ./2;
        xval = [ 1:(upperbound-lowerbound+1) ];
        linearmodel = grad_points .* xval + Datadiff(ii,lowerbound);
        %End copy code
        resnorm=zeros([numscans size(freqbounds,2)]);
        
        % To fit a Triple Gaussian

            % MM: Allowing peaks to vary individually seems to work better
            % than keeping the distance fixed (i.e., including J in the
            % function)
            
          % %x1-x6 for Glx, while x7-x9 for GABA --MGSaleh  
          % GaussModelInit = [  x(1)      x(2)  x(3)    x(4)    x(5)  x(6)    x(7)   x(8)  x(9)    x(10)      x(11)    x(12)  x(13)] --MGSaleh
            GaussModelInit = [maxinGABA  -400  3.725 maxinGABA -400  3.775 maxinGABA -90   3.02 -LinearInit constInit   0      0  ];
            lb = [0 -800 3.725-0.02 0 -800 3.775-0.02 0 -200 3.02-0.05 -40*maxinGABA -20000*maxinGABA -2000*maxinGABA -2000*maxinGABA];
            ub = [4000*maxinGABA -40 3.725+0.02 4000*maxinGABA -40 3.775+0.02 40000*maxinGABA -40 3.02+0.05 40*maxinGABA 1000*maxinGABA 1000*maxinGABA 1000*maxinGABA];

        
        options = optimset('lsqcurvefit');
        options = optimset(options,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',1e5);
        nlinopts = statset('nlinfit');
        nlinopts = statset(nlinopts, 'MaxIter', 1e5);
        
%                 figure(98)
%                 plot(freq(freqbounds),GABAGlxModel_area(GaussModelInit,freq(freqbounds)),freq(freqbounds),real(Datadiff(ii,freqbounds)))   
            [GaussModelParam(ii,:),resnorm,residg] = lsqcurvefit(@(xdummy,ydummy) GABAGlxModel_area(xdummy,ydummy), ...
                GaussModelInit, freq(freqbounds),real(Datadiff(ii,freqbounds)), ...
                lb,ub,options);
            residg = -residg;
            if(fit_method == FIT_NLINFIT)
                GaussModelInit = GaussModelParam(ii,:);
                
                
                for fit_iter = 1:1
                    [GaussModelParam(ii,:), residg, J, COVB, MSE] = nlinfit(freq(freqbounds), real(Datadiff(ii,freqbounds)), ... % J, COBV, MSE edited in
                        @(xdummy,ydummy) GABAGlxModel_area(xdummy,ydummy), ...
                        GaussModelInit, ...
                        nlinopts);
                    MRS_struct.out.(reg{kk}).fitparams_iter(fit_iter,:,ii) = GaussModelParam(ii,:);
                    GaussModelInit = GaussModelParam(ii,:);
                    ci = nlparci(GaussModelParam(ii,:), residg,'covar',COVB); %copied over
                end
                
%                 figure(99)
%                 plot(freq(freqbounds),GABAGlxModel_area(GaussModelInit,freq(freqbounds)),freq(freqbounds),real(Datadiff(ii,freqbounds)),freq(freqbounds),residg)
                
            end
            

            % Need to quantify GABA and Glx separately --MGSaleh

            % Glx fitting output -- MGSaleh
            MRS_struct.out.(reg{kk}).Glx.Area(ii) = (GaussModelParam(ii,1)./sqrt(-GaussModelParam(ii,2))*sqrt(pi)) + ...
                (GaussModelParam(ii,4)./sqrt(-GaussModelParam(ii,5))*sqrt(pi));
            Glxheight = max(GaussModelParam(ii,[1,4])); % The maximum Glx height is measured -- MGSaleh
            MRS_struct.out.(reg{kk}).Glx.FitError(ii)  =  100*std(residg([1:(midbound2-(lowerbound-1))],:))/Glxheight; % The residuals is from the whole range specified above for the the fitting --MGSaleh
            sigma = ((1/(2*(abs(GaussModelParam(ii,2))))).^(1/2)) + ((1/(2*(abs(GaussModelParam(ii,5))))).^(1/2));
            MRS_struct.out.(reg{kk}).Glx.FWHM(ii) =  abs( (2* MRS_struct.p.LarmorFreq) * sigma);
            MRS_struct.out.(reg{kk}).Glx.ModelFit(ii,:)=GaussModelParam(ii,:); % GaussModelParam(ii,[1:6]) are related to Glx -- MGSaleh;
            MRS_struct.out.(reg{kk}).Glx.resid(ii,:) = residg([1:(midbound2-(lowerbound-1))],:);
            MRS_struct.out.(reg{kk}).Glx.snr(ii) = Glxheight / std(residg([1:(midbound2-(lowerbound-1))],:));  % Added by MGSaleh            
            
      
            % GABA fitting output --MGSaleh
            MRS_struct.out.(reg{kk}).GABA.Area(ii) = (GaussModelParam(ii,7)./sqrt(-GaussModelParam(ii,8))*sqrt(pi));
            GABAheight = GaussModelParam(ii,7);
            MRS_struct.out.(reg{kk}).GABA.Height(ii) = GABAheight; 
            MRS_struct.out.(reg{kk}).GABA.FitError(ii)  =  100*std(residg([1:(upperbound-(midbound1-1))],:))/GABAheight; % The residuals is from the whole range specified above for the the fitting --MGSaleh
            sigma = ( 1 / (2 * (abs(GaussModelParam(ii,8)))) ).^(1/2);
            MRS_struct.out.(reg{kk}).GABA.FWHM(ii) =  abs( (2* MRS_struct.p.LarmorFreq) * sigma);
            MRS_struct.out.(reg{kk}).GABA.ModelFit(ii,:)=GaussModelParam(ii,:); % GaussModelParam(ii,7:9) are related to GABA  -- MGSaleh;
            MRS_struct.out.(reg{kk}).GABA.resid(ii,:) = residg([1:(upperbound-(midbound1-1))],:);
            MRS_struct.out.(reg{kk}).GABA.snr(ii) = GABAheight / std(residg([1:(upperbound-(midbound1-1))],:));  % Added by MGSaleh

            
    else
        error('Fitting MRS_struct.p.target not recognised');
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1A. Start up the output figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    fignum = 102;
    if(ishandle(fignum))
        close(fignum)
    end
    h=figure(fignum);
    set(h, 'Position', [100, 100, 1000, 707]);
    set(h,'Color',[1 1 1]);
    figTitle = ['GannetFit Output'];
    set(gcf,'Name',figTitle,'Tag',figTitle, 'NumberTitle','off');
    % GABA plot
    ha=subplot(2, 2, 1);
    % find peak of GABA plot... plot residuals above this...
    gabamin = min(real(Datadiff(ii,plotbounds)));
    gabamax = max(real(Datadiff(ii,plotbounds)));
    resmax = max(residg);
    residg = residg + gabamin - resmax;
    if strcmp(MRS_struct.p.target,'GABA')
        plot(freq(freqbounds),GaussModel_area(GaussModelParam(ii,:),freq(freqbounds)),'r',...
            freq(plotbounds),real(Datadiff(ii,plotbounds)), 'b', ...
            freq(freqbounds),residg,'k');
        set(gca,'XLim',[2.6 3.6]);
    elseif strcmp (target{trg},'GSH') %strcmp(MRS_struct.p.target,'GSH')
       freqrange=MRS_struct.spec.freq(freqbounds); 
       plot(MRS_struct.spec.freq, Datadiff(ii,:), 'b',freqrange, ...
            FiveGaussModel(FiveGaussModelParam(ii,:), freqrange),'r',freqrange, ...
            FiveGaussModel(GSHGaussModelParam(ii,:), freqrange),'g',...
            MRS_struct.spec.freq(freqbounds),residg,'k');
        set(gca,'XLim',[1.8 4.2]);
    elseif strcmp(MRS_struct.p.target,'Lac')
       freqrange=MRS_struct.spec.freq(freqbounds); 
       plot(MRS_struct.spec.freq, Datadiff(ii,:), 'b',freqrange, ...
            FourGaussModel(FourGaussModelParam(ii,:), freqrange),'r',freqrange, ...
            FourGaussModel(MMGaussModelParam(ii,:), freqrange),'r',...
            MRS_struct.spec.freq(freqbounds),residg,'k');
        set(gca,'XLim',[0.7 1.9]);
    elseif strcmp(MRS_struct.p.target,'Glx')
      plot(freq(freqbounds),DoubleGaussModel_area(GaussModelParam(ii,:),freq(freqbounds)),'r',...
                freq(plotbounds),real(Datadiff(ii,plotbounds)),'b', ...
                freq(freqbounds),residg,'k');
        set(gca,'XLim',[3.4 4.2]) 
    elseif strcmp (target{trg},'GABAGlx') %strcmp(MRS_struct.p.target,'GABAGlx')
      plot(freq(freqbounds),GABAGlxModel_area(GaussModelParam(ii,:),freq(freqbounds)),'r',...
                freq(plotbounds),real(Datadiff(ii,plotbounds)),'b', ...
                freq(freqbounds),residg,'k');
        set(gca,'XLim',[2.7 4.2])    
    end    
    if(strcmpi(MRS_struct.p.vendor,'Siemens'))
        legendtxt = regexprep(MRS_struct.gabafile{ii*2-1}, '_','-');
    else
        legendtxt = regexprep(MRS_struct.gabafile{ii}, '_','-');
    end
    title(legendtxt);
    set(gca,'XDir','reverse');
    
    if strcmp (target{trg},'GABA') 
        %%%%From here on is cosmetic - adding labels (and deciding where to).
        hgaba=text(3,gabamax/4,MRS_struct.p.target);
        set(hgaba,'horizontalAlignment', 'center');
        %determine values of GABA tail (below 2.8 ppm.
        z=abs(MRS_struct.spec.freq-2.79);%2.75
        upperbound=find(min(z)==z);
        tailtop=max(real(Datadiff(ii,upperbound:(upperbound+150))));
        tailbottom=min(real(Datadiff(ii,upperbound:(upperbound+150))));
        hgabares=text(2.8,min(residg),'residual');
        set(hgabares,'horizontalAlignment', 'left');
        text(2.8,tailtop+gabamax/20,'data','Color',[0 0 1]);
        text(2.8,tailbottom-gabamax/20,'model','Color',[1 0 0]);
    end
    
    
    if strcmp (target{trg},'GSH')
        %%%%From here on is cosmetic - adding labels (and deciding where to).
        hgaba=text(2.95,gabamax,target{trg});
        set(hgaba,'horizontalAlignment', 'center');
        %determine values of GABA tail (below 2.8 ppm.
        z=abs(MRS_struct.spec.freq-2.00);%2.75
        upperbound=find(min(z)==z);
        tailtop=max(real(Datadiff(ii,upperbound:(upperbound+150))));
        tailbottom=min(real(Datadiff(ii,upperbound:(upperbound+150))));
        hgabares=text(2.25,min(residg),'residual');
        set(hgabares,'horizontalAlignment', 'left');
        text(2.25,tailtop+gabamax/20,'data','Color',[0 0 1]);
        text(2.45,tailbottom-20*gabamax/20,'model','Color',[1 0 0]);
    end
    
    if strcmp(MRS_struct.p.target,'Glx')
        %%%%From here on is cosmetic - adding labels (and deciding where to).
        hgaba=text(3.8,gabamax/4,MRS_struct.p.target);
        set(hgaba,'horizontalAlignment', 'center');
        %determine values of GABA tail (below 2.8 ppm.
        z=abs(MRS_struct.spec.freq-2.79);%2.75
        upperbound=find(min(z)==z);
        tailtop=max(real(Datadiff(ii,upperbound:(upperbound+150))));
        tailbottom=min(real(Datadiff(ii,upperbound:(upperbound+150))));
        hgabares=text(3.5,min(residg),'residual');
        set(hgabares,'horizontalAlignment', 'left');
        text(3.5,tailtop+gabamax/20,'data','Color',[0 0 1]);
        text(3.5,tailbottom-gabamax/20,'model','Color',[1 0 0]);
    end
    
     if strcmp (target{trg},'GABAGlx')
        %%%%From here on is cosmetic - adding labels (and deciding where to).
        hgaba=text(3.2,gabamax/1.5,MRS_struct.p.target);
        set(hgaba,'horizontalAlignment', 'right');
        %determine values of GABA tail (below 2.8 ppm.
        z=abs(MRS_struct.spec.freq-2.79);%2.75
        upperbound=find(min(z)==z);
        tailtop=max(real(Datadiff(ii,upperbound:(upperbound+150))));
        tailbottom=min(real(Datadiff(ii,upperbound:(upperbound+150))));
        hgabares=text(2.78,min(residg),'residual');
        set(hgabares,'horizontalAlignment', 'left');
        text(2.78,tailtop+gabamax/20,'data','Color',[0 0 1]);
        text(2.78,tailbottom-gabamax/20,'model','Color',[1 0 0]);
    end
    set(gca,'YTick',[]);
    set(gca,'Box','off');
    set(gca,'YColor','white');
    
%% Water fitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2.  Water Fit 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(MRS_struct.p.Reference_compound,'H2O')
        T1=20;
        %estimate height and baseline from data
        [maxinWater, watermaxindex]=max(real(WaterData(ii,:)),[],2);
        waterbase = mean(real(WaterData(1:500))); % avg

        %Philips data do not phase well based on first point, so do a preliminary
        %fit, then adjust phase of WaterData accordingly     
        
        if(strcmpi(MRS_struct.p.vendor,'Philips'))
            %Run preliminary Fit of data
            LGModelInit = [maxinWater 20 freq(watermaxindex) 0.0 waterbase -50 ]; %works

            lblg = [0.01*maxinWater 1 4.6 0 0 -50 ];
            ublg = [40*maxinWater 100 4.8 0.000001 1 0 ];
            %Fit from 5.6 ppm to 3.8 ppm RE 110826
            z=abs(MRS_struct.spec.freq-5.6);
            waterlow=find(min(z)==z);
            z=abs(MRS_struct.spec.freq-3.8);
            waterhigh=find(min(z)==z);
            freqbounds=waterlow:waterhigh;
            % Do the water fit (Lorentz-Gauss)
            nlinopts = statset('nlinfit');
            nlinopts = statset(nlinopts, 'MaxIter', 1e5);
            [LGModelParam(ii,:),residw] = nlinfit(freq(freqbounds), real(WaterData(ii,freqbounds)),...
                @(xdummy,ydummy)	LorentzGaussModel(xdummy,ydummy),...
                LGModelInit, nlinopts);
            residw = -residw;
            
            %Then use this for phasing
            if ~MRS_struct.p.water_phase_correction % Do this only if the klose correction on water is off -- MGSaleh 2016 
            
            Eerror=zeros([120 1]);
            for jj=1:120
                Data=WaterData(ii,freqbounds)*exp(1i*pi/180*jj*3);
                Model=LorentzGaussModel(LGModelParam(ii,:),freq(freqbounds));
                Eerror(jj)=sum((real(Data)-Model).^2);
            end
            [number index]=min(Eerror);
            WaterData=WaterData*exp(1i*pi/180*index*3);
            
            end 
            
        end
    % x(1) = Amplitude of (scaled) Lorentzian
    % x(2) = 1 / hwhm of Lorentzian (hwhm = half width at half max)
    % x(3) = centre freq of Lorentzian
    % x(4) = linear baseline amplitude
    % x(5) = constant baseline amplitude
    % x(6) =  -1 / 2 * sigma^2  of gaussian
    LGModelInit = [maxinWater 20 4.7 0 waterbase -50 ]; %works
        lblg = [0.01*maxinWater 1 4.6 0 0 -50 ];
        ublg = [40*maxinWater 100 4.8 0.000001 1 0 ];
        %Fit from 5.6 ppm to 3.8 ppm RE 110826
        z=abs(MRS_struct.spec.freq-5.6);
        waterlow=find(min(z)==z);
        z=abs(MRS_struct.spec.freq-3.8);
        waterhigh=find(min(z)==z);
        freqbounds=waterlow:waterhigh;
        % Do the water fit (Lorentz-Gauss)
        % 111209 Always do the LSQCURV fitting - to initialise
            %Lorentz-Gauss Starters
            options = optimset('lsqcurvefit');
            options = optimset(options,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',10000);
            [LGModelParam(ii,:),residual(ii), residw] = lsqcurvefit(@(xdummy,ydummy) ...
                LorentzGaussModel(xdummy,ydummy),...
                LGModelInit, freq(freqbounds),real(WaterData(ii,freqbounds)),...
                lblg,ublg,options);
              residw = -residw;
            if(waterfit_method == FIT_NLINFIT)
                LGModelInit = LGModelParam(ii,:); % CJE 4 Jan 12   
                % nlinfit options
                nlinopts = statset('nlinfit');
                nlinopts = statset(nlinopts, 'MaxIter', 1e5);
                %This double fit doesn't seem to work too well with the GE
                %data... dig a little deeper
                LGPModelInit = [maxinWater 20 freq(watermaxindex) 0 waterbase -50 0];
                %figure(7)
                %plot(freq(freqbounds), real(WaterData(ii,freqbounds)),freq(freqbounds),LorentzGaussModelP(LGPModelInit,freq(freqbounds)))
                %figure(8)
                [LGPModelParam(ii,:),residw] = nlinfit(freq(freqbounds), real(WaterData(ii,freqbounds)),...
                    @(xdummy,ydummy)	LorentzGaussModelP(xdummy,ydummy),...
                    LGPModelInit, nlinopts);
                if(~strcmpi(MRS_struct.p.vendor,'GE')&&~strcmpi(MRS_struct.p.vendor,'Siemens'))
                    %remove phase and run again
                    WaterData(ii,:)=WaterData(ii,:)*exp(1i*LGPModelParam(ii,7));
                    LGPModelParam(ii,7)=0;
                    [LGPModelParam(ii,:),residw] = nlinfit(freq(freqbounds), real(WaterData(ii,freqbounds)),...
                    @(xdummy,ydummy)	LorentzGaussModelP(xdummy,ydummy),...
                    LGPModelParam(ii,:), nlinopts);
                end
                residw = -residw;
            end
        MRS_struct.out.(reg{kk}).Water.ModelParam(ii,:) = LGPModelParam(ii,:);

        hb=subplot(2, 2, 3);
        waterheight = LGPModelParam(ii,1);
        watmin = min(real(WaterData(ii,:)));
        watmax = max(real(WaterData(ii,:)));
        resmax = max(residw);
        MRS_struct.out.(reg{kk}).Water.FitError(ii)  =  100 * std(residw) / waterheight; %raee changed to residw
        residw = residw + watmin - resmax;
        stdevresidw=std(residw);
        
        % Measuring the FWHM of the water peak using the Model -- MGSaleh
        % 2016.
        MRS_struct.out.(reg{kk}).Water.FWHM(ii) = abs(fwhm(MRS_struct.spec.freq,real(LorentzGaussModelP(LGPModelParam(ii,:),MRS_struct.spec.freq))))*MRS_struct.p.LarmorFreq
       
        if strcmp (MRS_struct.p.target,'GABA') || strcmp (MRS_struct.p.target,'GABAGlx') % Added by MGSaleh 
        MRS_struct.out.(reg{kk}).GABA.IU_Error_w = (MRS_struct.out.(reg{kk}).GABA.FitError(ii) .^ 2 + ...
            MRS_struct.out.(reg{kk}).Water.FitError .^ 2 ) .^ 0.5;
        end
        
        if strcmp (MRS_struct.p.target,'Glx') || strcmp (MRS_struct.p.target,'GABAGlx')  % Added by MGSaleh
            MRS_struct.out.(reg{kk}).Glx.IU_Error_w = (MRS_struct.out.(reg{kk}).Glx.FitError(ii) .^ 2 + ...
            MRS_struct.out.(reg{kk}).Water.FitError .^ 2 ) .^ 0.5;
        end
        
        if strcmp ((target{trg}),'GSH')  % Added by MGSaleh
            MRS_struct.out.(reg{kk}).GSH.IU_Error_w = (MRS_struct.out.(reg{kk}).(target{trg}).FitError(ii) .^ 2 + ...
            MRS_struct.out.(reg{kk}).Water.FitError .^ 2 ) .^ 0.5;
        end
        
      
            
        plot(freq(freqbounds),real(LorentzGaussModelP(LGPModelParam(ii,:),freq(freqbounds))), 'r', ...
            freq(freqbounds),real(WaterData(ii,freqbounds)),'b', ...
            freq(freqbounds), residw, 'k');
        set(gca,'XDir','reverse');
        set(gca,'YTick',[]);
        set(gca,'Box','off');
        set(gca,'YColor','white');
        xlim([4.2 5.2]);
        %Add on some labels
        hwat=text(4.8,watmax/2,'Water');
        set(hwat,'horizontalAlignment', 'right');
        %Get the right vertical offset for the residual label
        z=abs(freq(freqbounds)-4.4);
        waterrlow=find(min(z)==z);
        z=abs(freq(freqbounds)-4.25);
        waterrhigh=find(min(z)==z);
        rlabelbounds=waterrlow:waterrhigh;
        labelfreq=freq(freqbounds);
        axis_bottom=axis;
        hwatres=text(4.4,max(min(residw(rlabelbounds))-0.05*watmax,axis_bottom(3)),'residual');
        set(hwatres,'horizontalAlignment', 'left');       
        %CJE fixes water baseline code - baseline model as before...
        WaterArea(ii)=sum(real(LorentzGaussModel(LGModelParam(ii,:),freq(freqbounds))) ...
      - BaselineModel(LGModelParam(ii,3:5),freq(freqbounds)),2);
        % convert watersum to integral
        MRS_struct.out.(reg{kk}).Water.Area(ii)=WaterArea(ii) * (freq(1) - freq(2));
        %MRS_struct.H20 = MRS_struct.out.WaterArea(ii) ./ std(residw); %This line doesn't make sense - commenting pending delete. RE
        %generate scaled spectrum (for plotting) CJE Jan2011
        Datadiff_scaled(ii,:) = Datadiff(ii,:) .* ...
            repmat((1 ./ MRS_struct.out.(reg{kk}).Water.Area(ii)), [1 32768]);
    
        % Concentration of GABA and Glx to water determined here. -- MGSaleh
        if strcmp (MRS_struct.p.target,'GABA') || strcmp (MRS_struct.p.target,'GABAGlx')
            [MRS_struct]=MRSGABAinstunits(MRS_struct, reg{kk},'GABA',ii);
        end
        
        if strcmp (MRS_struct.p.target,'Glx')|| strcmp (MRS_struct.p.target,'GABAGlx')
            [MRS_struct]=MRSGABAinstunits(MRS_struct,reg{kk},'Glx', ii);
        end
       
        if strcmp ((target{trg}),'GSH')
            [MRS_struct]=MRSGABAinstunits(MRS_struct,reg{kk},(target{trg}), ii);
        end

            
    end
    
    
    
%% Creatine fitting
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3.  Cr Fit 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
Cr_OFF=Dataoff(ii,:);        
%Fit CHo and Cr
    ChoCrFitLimLow=2.6;
    ChoCrFitLimHigh=3.6;           
    %Still need ranges for Creatine align plot
    z=abs(MRS_struct.spec.freq-ChoCrFitLimHigh);
    cclb=find(min(z)==z);
    z=abs(MRS_struct.spec.freq-ChoCrFitLimLow);
    ccub=find(min(z)==z);
    freqrangecc=MRS_struct.spec.freq(cclb:ccub);
    %Do some detective work to figure out the initial parameters
    ChoCrMeanSpec = Cr_OFF(cclb:ccub).';
    Baseline_offset=real(ChoCrMeanSpec(1)+ChoCrMeanSpec(end))/2;
    Width_estimate=0.05;%ppm
    Area_estimate=(max(real(ChoCrMeanSpec))-min(real(ChoCrMeanSpec)))*Width_estimate*4;
    ChoCr_initx = [ Area_estimate Width_estimate 3.02 0 Baseline_offset 0 1].*[1 (2*MRS_struct.p.LarmorFreq) MRS_struct.p.LarmorFreq (180/pi) 1 1 1];     
    ChoCrMeanSpecFit(ii,:) = FitChoCr(freqrangecc, ChoCrMeanSpec, ChoCr_initx,MRS_struct.p.LarmorFreq);
    MRS_struct.out.ChoCrMeanSpecFit(ii,:) = ChoCrMeanSpecFit(ii,:)./[1 (2*MRS_struct.p.LarmorFreq) MRS_struct.p.LarmorFreq (180/pi) 1 1 1];

        %Initialise fitting pars
        z=abs(MRS_struct.spec.freq-3.12);
        lb=find(min(z)==z);
        z=abs(MRS_struct.spec.freq-2.72);
        ub=find(min(z)==z);
        Cr_initx = [max(real(Cr_OFF(lb:ub))) 0.05 3.0 0 0 0 ];
        freqrange = MRS_struct.spec.freq(lb:ub);
        %Then use the same function as the Cr Fit in GannetLoad
        nlinopts=statset('nlinfit');
        nlinopts = statset(nlinopts, 'MaxIter', 1e5, 'Display','Off');
        [CrFitParams(ii,:), residCr] = nlinfit(freqrange, real(Cr_OFF(lb:ub)), ...
           @(xdummy, ydummy) LorentzModel(xdummy, ydummy),Cr_initx, nlinopts);
        Crheight = CrFitParams(ii,1);
        Crmin = min(real(Cr_OFF(lb:ub)));
        Crmax = max(real(Cr_OFF(lb:ub)));
        resmaxCr = max(residCr);
        stdresidCr = std(residCr);
        MRS_struct.out.Creatine.FitError(ii)  =  100 * stdresidCr / Crheight;
       
        if strcmp (MRS_struct.p.target,'GABA') || strcmp (MRS_struct.p.target,'GABAGlx') % Added by MGSaleh
        MRS_struct.out.(reg{kk}).GABA.IU_Error_cr(ii) = (MRS_struct.out.(reg{kk}).GABA.FitError(ii) .^ 2 + ...
            MRS_struct.out.Creatine.FitError(ii) .^ 2 ) .^ 0.5;
        end
        
        if strcmp (MRS_struct.p.target,'Glx') || strcmp (MRS_struct.p.target,'GABAGlx') % Added by MGSaleh
             MRS_struct.out.(reg{kk}).Glx.IU_Error_cr(ii) = (MRS_struct.out.(reg{kk}).Glx.FitError(ii) .^ 2 + ...
            MRS_struct.out.Creatine.FitError(ii) .^ 2 ) .^ 0.5;
        end
        
        if strcmp ((target{trg}),'GSH')  % Added by MGSaleh 2017
             MRS_struct.out.(reg{kk}).(target{trg}).IU_Error_cr(ii) = (MRS_struct.out.(reg{kk}).(target{trg}).FitError(ii) .^ 2 + ...
            MRS_struct.out.Creatine.FitError(ii) .^ 2 ) .^ 0.5;
        end
            
        %MRS_struct.out.CrArea(ii)=sum(real(LorentzModel(CrFitParams(ii,:),freqrange)-LorentzModel([0 CrFitParams(ii,2:end)],freqrange))) * (freq(1) - freq(2));
        
        MRS_struct.out.(reg{kk}).Creatine.Area(ii)=sum(real(TwoLorentzModel([MRS_struct.out.ChoCrMeanSpecFit(ii,1:(end-1)) 0],freqrangecc)-TwoLorentzModel([0 MRS_struct.out.ChoCrMeanSpecFit(ii,2:(end-1)) 0],freqrangecc))) * (freq(1) - freq(2));
        MRS_struct.out.(reg{kk}).Choline.Area(ii)=sum(real(TwoLorentzModel([MRS_struct.out.ChoCrMeanSpecFit(ii,1:(end))],freqrangecc)-TwoLorentzModel([MRS_struct.out.ChoCrMeanSpecFit(ii,1:(end-1)) 0],freqrangecc))) * (freq(1) - freq(2));
        MRS_struct.out.CrFWHMHz(ii)= ChoCrMeanSpecFit(ii,2);
        
        if strcmp (MRS_struct.p.target,'GABA') || strcmp (MRS_struct.p.target,'GABAGlx') % Added by MGSaleh
            MRS_struct.out.(reg{kk}).GABA.GABAconcCr(ii)=MRS_struct.out.(reg{kk}).GABA.Area(ii)./MRS_struct.out.(reg{kk}).Creatine.Area(ii);   
            MRS_struct.out.(reg{kk}).GABA.GABAconcCho(ii)=MRS_struct.out.(reg{kk}).GABA.Area(ii)./MRS_struct.out.(reg{kk}).Choline.Area(ii);    
        end
        
        if strcmp (MRS_struct.p.target,'Glx') || strcmp (MRS_struct.p.target,'GABAGlx') % Added by MGSaleh
            MRS_struct.out.(reg{kk}).Glx.GlxconcCr(ii)=MRS_struct.out.(reg{kk}).Glx.Area(ii)./MRS_struct.out.(reg{kk}).Creatine.Area(ii);   
            MRS_struct.out.(reg{kk}).Glx.GlxconcCho(ii)=MRS_struct.out.(reg{kk}).Glx.Area(ii)./MRS_struct.out.(reg{kk}).Choline.Area(ii); 
            
        end
        
        if strcmp ((target{trg}),'GSH') % Added by MGSaleh
            MRS_struct.out.(reg{kk}).(target{trg}).GSHconcCr(ii)=MRS_struct.out.(reg{kk}).(target{trg}).Area(ii)./MRS_struct.out.(reg{kk}).Creatine.Area(ii);   
            MRS_struct.out.(reg{kk}).(target{trg}).GSHconcCho(ii)=MRS_struct.out.(reg{kk}).(target{trg}).Area(ii)./MRS_struct.out.(reg{kk}).Choline.Area(ii);    
        end
        
%% Output parameters and display options        
                
        %alter resid Cr for plotting.
        residCr = residCr + Crmin - resmaxCr;
        if strcmp(MRS_struct.p.Reference_compound,'H2O')
            %Plot the Cr fit
            h2=subplot(2, 2, 4);
            %debugging changes
            plot(freqrangecc,real(TwoLorentzModel(MRS_struct.out.ChoCrMeanSpecFit(ii,:),freqrangecc)), 'r', ...
                freqrangecc,real(TwoLorentzModel([MRS_struct.out.ChoCrMeanSpecFit(ii,1:(end-1)) 0],freqrangecc)), 'r', ...
                MRS_struct.spec.freq,real(Cr_OFF(:)),'b', ...
                freqrange, residCr, 'k');
            set(gca,'XDir','reverse');
            set(gca,'YTick',[],'Box','off');
            xlim([2.6 3.6]);
            set(gca,'YColor','white');
            hcr=text(2.94,Crmax*0.75,'Creatine');
            set(hcr,'horizontalAlignment', 'left')
            %Transfer Cr plot into insert
            subplot(2,2,3)
            [h_m h_i]=inset(hb,h2);
            set(h_i,'fontsize',6);
            insert=get(h_i,'pos');
            axi=get(hb,'pos');
            set(h_i,'pos',[axi(1)+axi(3)-insert(3) insert(2:4)]);
            %Add labels
            hwat=text(4.8,watmax/2,'Water');
            set(hwat,'horizontalAlignment', 'right')
            set(h_m,'YTickLabel',[]);
            set(h_m,'XTickLabel',[]);
            set(gca,'Box','off')
            set(gca,'YColor','white');
        else
            %Plot the Cr fit
            hb=subplot(2, 2, 3);
            %debugging changes
            %plot(freqrange,real(LorentzModel(CrFitParams(ii,:),freqrange)), 'r', ...
            %    MRS_struct.spec.freq,real(Cr_OFF(:)),'b', ...
            %    freqrange, residCr, 'k');            
            plot(freqrangecc,real(TwoLorentzModel(MRS_struct.out.ChoCrMeanSpecFit(ii,:),freqrangecc)), 'r', ...
                freqrangecc,real(TwoLorentzModel([MRS_struct.out.ChoCrMeanSpecFit(ii,1:(end-1)) 0],freqrangecc)), 'r', ...
                MRS_struct.spec.freq,real(Cr_OFF(:)),'b', ...
                freqrange, residCr, 'k');
            
            set(gca,'XDir','reverse');
            set(gca,'YTick',[]);
            xlim([2.6 3.6]);
            z=abs(freq(lb:ub)-3.12);
            crlow=find(min(z)==z);
            z=abs(freq(lb:ub)-2.9);
            crhigh=find(min(z)==z);
            crlabelbounds=crlow:crhigh;
            hcres=text(3.12,max(residCr(crlabelbounds))+0.05*Crmax,'residual');
            set(hcres,'horizontalAlignment', 'left');
            hcdata=text(2.8,0.3*Crmax,'data','Color',[0 0 1]);
            hcmodel=text(2.8,0.2*Crmax,'model','Color',[1 0 0]);
            text(2.94,Crmax*0.75,'Creatine');
        end

    % GABA fitting information
    if(strcmp(MRS_struct.p.AlignTo,'no')~=1)
        tmp2 = '1';
    else
        tmp2 = '0';
    end
    if fit_method == FIT_NLINFIT
        tmp3 = 'NLINFIT, ';
    else
        tmp3 = 'LSQCURVEFIT, ';
    end
    if waterfit_method == FIT_NLINFIT
        tmp4 = [tmp3 'NLINFIT'];
    else
        tmp4 = [tmp3 'LSQCURVEFIT' ];
    end


    %and running the plot
    if any(strcmp('mask',fieldnames(MRS_struct))) == 1
    h=subplot(2,2,2);
    p = get(h,'pos'); % get position of axes
    set(h,'pos',[0.52 0.52 0.42 0.42]) % move the axes slightly

    input=MRS_struct.mask.img(MRS_struct.ii,:,1:round(size(MRS_struct.mask.img,3)/3));
   
    imagesc(squeeze(MRS_struct.mask.img(MRS_struct.ii,:,1:round(size(MRS_struct.mask.img,3)/3))));
    colormap('gray');
    caxis([0 1])
    axis equal;
    axis tight;
    axis off;
    subplot(2,2,4,'replace')
    else
    subplot(2,2,2)        
    end
    
    axis off
      if strcmp(MRS_struct.p.vendor,'Siemens')
         tmp = [ 'filename    : ' MRS_struct.gabafile{ii*2-1} ];
     else
        tmp = [ 'filename    : ' MRS_struct.gabafile{ii} ];
     end
    tmp = regexprep(tmp, '_','-');
    
    var_pos=0.9; % A variable added to determine position of the printout on the PDF/matlab-output window -- Added by MGSaleh
    
    text(0,var_pos, tmp, 'FontName', 'Helvetica');
    if isfield(MRS_struct.p,'voxsize')
    tmp =       [ num2str(MRS_struct.p.Navg(ii)) ' averages of a ' num2str(MRS_struct.p.voxsize(ii,1)*MRS_struct.p.voxsize(ii,2)*MRS_struct.p.voxsize(ii,3)*.001) ' ml voxel'];
    else
    tmp =       [ num2str(MRS_struct.p.Navg(ii)) ' averages'];
    end
    text(0,var_pos-0.1, tmp, 'FontName', 'Helvetica');
    %Remove this - more useful to add in Cr fWHM at a later date
    %tmp = sprintf('GABA+ FWHM   : %.2f Hz', MRS_struct.out.GABAFWHM(ii) );
    %text(0,0.7, tmp);
    if isfield(MRS_struct.p,'voxsize')
             SNRfactor=round(sqrt(MRS_struct.p.Navg(ii))*MRS_struct.p.voxsize(ii,1)*MRS_struct.p.voxsize(ii,2)*MRS_struct.p.voxsize(ii,3)*.001);
             tmp = ['SNR factor         :  '  num2str(SNRfactor)];
             text(0,var_pos-0.2, tmp, 'FontName', 'Helvetica');
    end
    
    % Some changes to the printout to accomodate GABAGlx fitting output
    if strcmp (MRS_struct.p.target, 'GABA') % Added by MGSaleh
        tmp = sprintf('GABA+ Area : %.3g', MRS_struct.out.(reg{kk}).GABA.Area(ii));
    elseif strcmp (MRS_struct.p.target, 'Glx') % Added by MGSaleh
        tmp = sprintf('Glx Area : %.3g', MRS_struct.out.(reg{kk}).Glx.Area(ii));
    elseif strcmp ((target{trg}), 'GABAGlx') % Added by MGSaleh
        tmp = sprintf('GABA+/Glx Areas : %.3g/%.3g', MRS_struct.out.(reg{kk}).GABA.Area(ii),MRS_struct.out.(reg{kk}).Glx.Area(ii));

    end
    
    if strcmp ((target{trg}), 'GSH')  % Added by MGSaleh 2017
        tmp = sprintf('GSH Area : %.3g', MRS_struct.out.(reg{kk}).(target{trg}).Area(ii));
    end
    text(0,var_pos-0.3, tmp);

    % Made some changes to FitErr printout to accommodate GABAGlx printouts
    if strcmp(MRS_struct.p.Reference_compound,'H2O')
        tmp = sprintf('FWHM of Water/Cr: %.1f/%.1f Hz ', MRS_struct.out.(reg{kk}).Water.FWHM(ii),MRS_struct.out.CrFWHMHz(ii)  );
        text(0,var_pos-0.5, tmp, 'FontName', 'Helvetica');
        tmp = sprintf('H_2O/Cr Area : %.3g/%.3g ', MRS_struct.out.(reg{kk}).Water.Area(ii),MRS_struct.out.(reg{kk}).Creatine.Area(ii) );
        text(0,var_pos-0.4, tmp, 'FontName', 'Helvetica');
        
        if strcmp (MRS_struct.p.target, 'GABA') % Added by MGSaleh
            tmp = sprintf('%.1f, %.1f ',  MRS_struct.out.GABA.IU_Error_w(ii),  MRS_struct.out.GABA.IU_Error_cr(ii));
            tmp = [tmp '%'];
            tmp = ['FitErr (H/Cr)   : ' tmp];
        elseif strcmp (MRS_struct.p.target, 'Glx') % Added by MGSaleh
            tmp = sprintf('%.1f, %.1f ',  MRS_struct.out.Glx.IU_Error_w(ii),  MRS_struct.out.Glx.IU_Error_cr(ii));
            tmp = [tmp '%'];
            tmp = ['FitErr (H/Cr)   : ' tmp];
        elseif strcmp ((target{trg}), 'GABAGlx') % Added by MGSaleh
            tmp = sprintf('GABA(H/Cr):   %.1f/%.1f;   Glx(H/Cr):   %.1f/%.1f ',  MRS_struct.out.(reg{kk}).GABA.IU_Error_w(ii),  MRS_struct.out.(reg{kk}).GABA.IU_Error_cr(ii), MRS_struct.out.(reg{kk}).Glx.IU_Error_w(ii),  MRS_struct.out.(reg{kk}).Glx.IU_Error_cr(ii));
            tmp = ['FitErr(%) --  ' tmp];
        end
        
        if strcmp ((target{trg}), 'GSH')  % Added by MGSaleh 2017
            tmp = sprintf('GSH(H/Cr):   %.1f/%.1f',  MRS_struct.out.(reg{kk}).(target{trg}).IU_Error_w(ii),  MRS_struct.out.(reg{kk}).(target{trg}).IU_Error_cr(ii));
            tmp = ['FitErr(%) --  ' tmp];
        end

               
        if strcmp (MRS_struct.p.target, 'GABA') % Added by MGSaleh
            text(0,var_pos-0.6, tmp, 'FontName', 'Helvetica');
            tmp = [MRS_struct.p.target sprintf( '+/H_2O: %.3f inst. units.', MRS_struct.out.(reg{kk}).GABA.conciu(ii) )];
            text(0,var_pos-0.7, tmp, 'FontName', 'Helvetica');
            tmp = [ MRS_struct.p.target sprintf('+/Cr i.r.: %.3f', MRS_struct.out.(reg{kk}).GABA.GABAconcCr(ii) )]; 
        
        elseif strcmp (MRS_struct.p.target, 'Glx') % Added by MGSaleh
            text(0,var_pos-0.6, tmp, 'FontName', 'Helvetica');
            tmp = [MRS_struct.p.target sprintf( '/H_2O: %.3f inst. units.', MRS_struct.out.(reg{kk}).Glx.conciu(ii) )];
            text(0,var_pos-0.7, tmp, 'FontName', 'Helvetica');
            tmp = [ MRS_struct.p.target sprintf('/Cr i.r.: %.3f', MRS_struct.out.(reg{kk}).Glx.GlxconcCr(ii) )];  
            
        elseif strcmp ((target{trg}), 'GABAGlx') % Added by MGSaleh
            text(0,var_pos-0.6, tmp, 'FontName', 'Helvetica');
            tmp = [ sprintf( '[GABA+]/H_2O:   %.3f;   [Glx]/H_2O:   %.3f inst. units.', MRS_struct.out.(reg{kk}).GABA.conciu(ii), MRS_struct.out.(reg{kk}).Glx.conciu(ii) )];
            text(0,var_pos-0.7, tmp, 'FontName', 'Helvetica');
            tmp = [ sprintf( '[GABA+]/Cr i.r.:   %.3f;   [Glx]/Cr i.r.:   %.3f', MRS_struct.out.(reg{kk}).GABA.GABAconcCr(ii),MRS_struct.out.(reg{kk}).Glx.GlxconcCr(ii) )];
            
        end
        
        if strcmp ((target{trg}), 'GSH')  % Added by MGSaleh 2017
            text(0,var_pos-0.6, tmp, 'FontName', 'Helvetica');
            tmp = [ sprintf( '[GSH]/H_2O:   %.3f inst. units.', MRS_struct.out.(reg{kk}).(target{trg}).conciu(ii) )];
            text(0,var_pos-0.7, tmp, 'FontName', 'Helvetica');
            tmp = [ sprintf( '[GSH]/Cr i.r.:   %.3f', MRS_struct.out.(reg{kk}).(target{trg}).GSHconcCr(ii) )];
            
        end
        
        text(0,var_pos-0.8, tmp, 'FontName', 'Helvetica');
        tmp =       [ 'Ver(Load/Fit): ' MRS_struct.versionload  ',' MRS_struct.versionfit];
        text(0,var_pos-0.9, tmp, 'FontName', 'Helvetica');
        %tmp =        [MRS_struct.p.target ', Water fit alg. :' tmp4 ];
        %text(0,-0.1, tmp, 'FontName', 'Helvetica');
    else
       
        tmp = sprintf('Cr Area      : %.4f', MRS_struct.out.(reg{kk}).Creatine.Area(ii) );
        text(0,var_pos-0.4, tmp, 'FontName', 'Helvetica');
        
        if strcmp ((target{trg}),'GABA')       
            tmp = sprintf('GABA(/Cr): %.1f',  MRS_struct.out.(reg{kk}).GABA.IU_Error_cr(ii));
            tmp = ['FitErr(%) -- ' tmp ];
            text(0,var_pos-0.5, tmp, 'FontName', 'Helvetica');
            
        elseif strcmp ((target{trg}), 'Glx')
            tmp = sprintf('Glx(/Cr) : %.1f',  MRS_struct.out.(reg{kk}).Glx.IU_Error_cr(ii));
            tmp = ['FitErr(%) -- ' tmp ];
            tmp = ['FitErr(%) -- ' tmp ];
            text(0,var_pos-0.5, tmp, 'FontName', 'Helvetica');
            
        elseif strcmp ((target{trg}), 'GABAGlx') % Added by MGSaleh
            tmp = sprintf('GABA(/Cr):   %.1f;   Glx(/Cr):  %.1f ',  MRS_struct.out.(reg{kk}).GABA.IU_Error_cr(ii), MRS_struct.out.(reg{kk}).Glx.IU_Error_cr(ii));
            tmp = ['FitErr(%) -- ' tmp ];
            text(0,var_pos-0.5, tmp, 'FontName', 'Helvetica');
        
        else strcmp ((target{trg}), 'GSH') % Added by MGSaleh
            tmp = sprintf('GSH(/Cr) : %.1f',  MRS_struct.out.(reg{kk}).GSH.IU_Error_cr(ii));
            tmp = ['FitErr(%) -- ' tmp ];
            text(0,var_pos-0.5, tmp, 'FontName', 'Helvetica');
        end
        
        
        if strcmp ((target{trg}),'GABA')       
            tmp = [MRS_struct.p.target sprintf( '+/Cr i.r.: %.4f', MRS_struct.out.(reg{kk}).GABA.GABAconcCr(ii) )];
            text(0,var_pos-0.6, tmp, 'FontName', 'Helvetica');
        elseif strcmp ((target{trg}), 'Glx')
            tmp = [MRS_struct.p.target sprintf( '/Cr i.r.: %.4f', MRS_struct.out.(reg{kk}).Glx.GlxconcCr(ii) )];
            text(0,var_pos-0.6, tmp, 'FontName', 'Helvetica');
        elseif strcmp ((target{trg}), 'GABAGlx') % Added by MGSaleh
            tmp = [sprintf( '[GABA+]/Cr i.r.:   %.3f;   [Glx]/Cr i.r.:   %.3f', MRS_struct.out.(reg{kk}).GABA.GABAconcCr(ii),MRS_struct.out.(reg{kk}).Glx.GlxconcCr(ii) )];
            text(0,var_pos-0.6, tmp, 'FontName', 'Helvetica');
        end
        
        if strcmp ((target{trg}), 'GSH') % Added by MGSaleh 2017
           tmp = [MRS_struct.p.target sprintf( '/Cr i.r.: %.4f', MRS_struct.out.(reg{kk}).(target{trg}).GSHconcCr(ii) )];
           text(0,var_pos-0.6, tmp, 'FontName', 'Helvetica');
        end
        
        
        tmp =       [ 'Ver(Load/Fit): ' MRS_struct.versionload ',' MRS_struct.versionfit];
        text(0,var_pos-.7, tmp, 'FontName', 'Helvetica');
        %tmp =        [MRS_struct.p.target ', Water fit alg. :' tmp4 ];
        %text(0,0.0, tmp, 'FontName', 'Helvetica');
    end
    %Add Gannet logo
    if any(strcmp('mask',fieldnames(MRS_struct))) == 1
    
    subplot(2,2,4)
    else
    subplot(2,2,4,'replace')
    end
    axis off;
    script_path=which('GannetFit');
    Gannet_circle_white=[script_path(1:(end-12)) '/GANNET_circle_white.jpg'];
    A_2=imread(Gannet_circle_white);
    hax=axes('Position',[0.80, 0.05, 0.15, 0.15]);
    image(A_2);axis off; axis square;

    %%%%  Save EPS %%%%%
    if strcmp(MRS_struct.p.vendor,'Siemens')
    pfil_nopath = MRS_struct.gabafile{ii*2-1};
    else
    pfil_nopath = MRS_struct.gabafile{ii};
    end
    %for philips .data
    if(strcmpi(MRS_struct.p.vendor,'Philips_data'))
        fullpath = MRS_struct.gabafile{ii};
        %         fullpath = regexprep(fullpath, '\./', ''); NP edit out.
        %         see below
        %         fullpath = regexprep(fullpath, '/', '_');
        fullpath = regexprep(fullpath, '.data', '_data');
        fullpath = regexprep(fullpath, '\', '_');
        fullpath = regexprep(fullpath, '/', '_');
        %NP edit 02012013
        %Previous code somehow didn't run when running from hierarchical
        %folder (e.g. GABA_file = '.\name\MRI\raw.data) I got an error when Gannet tried to save the pdf for
        %.data file. E.g. ??? Error using ==> saveas at 115 Invalid or missing path: ./MRSfit_140102/.\7011-0124\MRI\raw_008.data.pdf
        %So it obviously didn't rewrite the path properly for the pdf here, but it IS important to get both folder and filename
        %as a lot of the .data files have similar names (e.g.
        %%raw_001.data). This change works for me for now, might not
        %%be most elegant
        
    end
    tmp = strfind(pfil_nopath,'/');
    tmp2 = strfind(pfil_nopath,'\');
    if(tmp)
        lastslash=tmp(end);
    elseif (tmp2)
        %maybe it's Windows...
        lastslash=tmp2(end);
    else
        % it's in the current dir...
        lastslash=0;
    end
    if(strcmpi(MRS_struct.p.vendor,'Philips'))
        tmp = strfind(pfil_nopath, '.sdat');
        tmp1= strfind(pfil_nopath, '.SDAT');
        if size(tmp,1)>size(tmp1,1)
            dot7 = tmp(end); % just in case there's another .sdat somewhere else...
        else
            dot7 = tmp1(end); % just in case there's another .sdat somewhere else...
        end
    elseif(strcmpi(MRS_struct.p.vendor,'GE'))
        tmp = strfind(pfil_nopath, '.7');
        dot7 = tmp(end); % just in case there's another .7 somewhere else...
    elseif(strcmpi(MRS_struct.p.vendor,'Philips_data'))
        tmp = strfind(pfil_nopath, '.data');
        dot7 = tmp(end); % just in case there's another .data somewhere else...
    elseif(strcmpi(MRS_struct.p.vendor,'Siemens'))
        tmp = strfind(pfil_nopath, '.rda');
        dot7 = tmp(end); % just in case there's another .data somewhere else...
    elseif(strcmpi(MRS_struct.p.vendor,'Siemens_twix'))
        tmp = strfind(pfil_nopath, '.dat');
        dot7 = tmp(end); % just in case there's another .dat somewhere else...
    end
    pfil_nopath = pfil_nopath( (lastslash+1) : (dot7-1) );
    if sum(strcmp(listfonts,'Helvetica'))>0
           set(findall(h,'type','text'),'FontName','Helvetica')
           set(ha,'FontName','Helvetica')
           set(hb,'FontName','Helvetica')
    end
    
%% Save pdf output

    set(gcf, 'PaperUnits', 'inches');
    set(gcf,'PaperSize',[11 8.5]);
    set(gcf,'PaperPosition',[0 0 11 8.5]);
    if(strcmpi(MRS_struct.p.vendor,'Philips_data'))
        pdfname=[ epsdirname '/' fullpath '-' target{trg} '.pdf' ];
    else
        pdfname=[ epsdirname '/' pfil_nopath  '-' target{trg} '.pdf' ];
    end
    %epsdirname
    if(exist(epsdirname,'dir') ~= 7)
        epsdirname
        mkdir(epsdirname)
    end
    saveas(gcf, pdfname);
    if(ii==numscans)
    if((MRS_struct.p.mat) == 1)
       if(strcmpi(MRS_struct.p.vendor,'Philips_data'))
       matname=[ epsdirname '/' 'MRS_struct'  '.mat' ];
       else
        matname =[ epsdirname '/' 'MRS_struct' '.mat' ];
       end
       
         save(matname,'MRS_struct'); 
    end
    end
    
    
%140116: ADH reorder structure


      if(isfield(MRS_struct, 'mask') == 1)

       if(isfield(MRS_struct, 'waterfile') == 1)
            structorder = {'versionload', 'versionfit', 'ii', ...
                'gabafile', 'waterfile', 'p', 'fids', 'spec', 'out', 'mask'};
        else 
             structorder = {'versionload', 'versionfit','ii', ...
                 'gabafile', 'p', 'fids', 'spec', 'out', 'mask'};
       end
       
      else
             
       if(isfield(MRS_struct, 'waterfile') == 1)
            structorder = {'versionload', 'versionfit', 'ii', ...
                'gabafile', 'waterfile', 'p', 'fids', 'spec', 'out'};
        else 
             structorder = {'versionload', 'versionfit','ii', ...
                 'gabafile', 'p', 'fids', 'spec', 'out'};
       end
      end
        


    
% Dec 09: based on FitSeries.m:  Richard's GABA Fitting routine
%     Fits using GaussModel
% Feb 10: Change the quantification method for water.  Regions of poor homogeneity (e.g. limbic)
%     can produce highly asymetric lineshapes, which are fitted poorly.  Don't fit - integrate
%     the water peak.
% March 10: 100301
%           use MRS_struct to pass loaded data data, call MRSGABAinstunits from here.
%           scaling of fitting to sort out differences between original (RE) and my analysis of FEF data
%           change tolerance on gaba fit
% 110308:   Keep definitions of fit functions in MRSGABAfit, rather
%               than in separate .m files
%           Ditto institutional units calc
%           Include FIXED version of Lorentzian fitting
%           Get Navg from struct (need version 110303, or later of
%               MRSLoadPfiles
%           rejig the output plots - one fig per scan.
% 110624:   set parmeter to choose fitting routine... for awkward spectra
%           report fit error (100*stdev(resid)/gabaheight), rather than "SNR"
%           can estimate this from confidence interval for nlinfit - need
%               GABA and water estimates

% 111111:   RAEE To integrate in Philips data, which doesn't always have
% water spectr, we need to add in referenceing to Cr... through
% MRS_struct.p.Reference_compound
% 140115: MRS_struct.p.Reference_compound is now 
%   MRS_struct.p.Reference compound
% 
%111214 integrating CJE's changes on water fitting (pre-init and revert to
%linear bseline). Also investigating Navg(ii)

        end 
end

end          % end of MRSGABAfit



%%%%%%%%%%%%%%%%%%%%%%%% GAUSS MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = GaussModel_area(x,freq)
%% Function for Gauss Model 


% x(1) = gaussian amplitude
% x(2) = 1/(2*sigma^2)
% x(3) = centre freq of peak
% x(4) = amplitude of linear baseline
% x(5) = constant amplitude offset

%F = x(1)*sqrt(-x(2)/pi)*exp(x(2)*(freq-x(3)).*(freq-x(3)))+x(4)*(freq-x(3))+x(5);
F = x(1)*exp(x(2)*(freq-x(3)).*(freq-x(3)))+x(4)*(freq-x(3))+x(5);



%%%%%%%%%%%%%%%%  OLD LORENTZGAUSSMODEL %%%%%%%%%%%%%%%%%%%%
%function F = LorentzGaussModel(x,freq)
%Lorentzian Model multiplied by a Gaussian.  gaussian width determined by
%x(6). x(7) determines phase.
%F = ((ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)*x(1))*cos(x(7))+(ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)*x(2).*(freq-x(3)))*sin(x(7))).*(exp(x(6)*(freq-x(3)).*(freq-x(3))))+x(4)*(freq-x(3))+x(5);



%%%%%%%%%%%%%%%%  LORENTZGAUSSMODEL %%%%%%%%%%%%%%%%%%%%
function F = LorentzGaussModel(x,freq)
%% Function for LorentzGaussModel Model 


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

%F =((ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)*x(1))*cos(x(7))+(ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)*x(2).*(freq-x(3)))*sin(x(7))).*(exp(x(6)*(freq-x(3)).*(freq-x(3))))+x(4)*(freq-x(3))+x(5);
% remove phasing
F = (x(1)*ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1))  ...
    .* (exp(x(6)*(freq-x(3)).*(freq-x(3)))) ... % gaussian
    + x(4)*(freq-x(3)) ... % linear baseline
    +x(5); % constant baseline

%%%%%%%%%%%%%%%%  NEW LORENTZGAUSSMODEL WITH PHASE%%%%%%%%%%%%%%%%%%%%
function F = LorentzGaussModelP(x,freq)
%% Function for LorentzGaussModel Model with Phase 


% CJE 24Nov10 - removed phase term from fit - this is now dealt with
% by the phasing of the water ref scans in MRSLoadPfiles
%Lorentzian Model multiplied by a Gaussian.
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

%F =((ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)*x(1))*cos(x(7))+(ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)*x(2).*(freq-x(3)))*sin(x(7))).*(exp(x(6)*(freq-x(3)).*(freq-x(3))))+x(4)*(freq-x(3))+x(5);
% remove phasing
F = ((cos(x(7))*x(1)*ones(size(freq))+sin(x(7)*x(1)*x(2)*(freq-x(3))))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1))  ...
    .* (exp(x(6)*(freq-x(3)).*(freq-x(3)))) ... % gaussian
    + x(4)*(freq-x(3)) ... % linear baseline
    +x(5); % constant baseline

  
%%%%%%%%%%%%%%%%%%%%%%%% DOUBLE GAUSS MODEL (MM: 150211) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = DoubleGaussModel_area(x,freq)
%% Function for DoubleGaussModel Model  


% A = 5/127.7; % for scalar coupling for Glu-C2 multiplet at 3.74 ppm (5 Hz seems to work best for in vivo data)
% F = x(1)*exp(x(2)*(freq-x(3)+A).*(freq-x(3)+A)) ...
%     +x(1)*exp(x(2)*(freq-x(3)-A).*(freq-x(3)-A)) ...
%     +x(4)*(freq-x(3))+x(5);

% Two Gaussians
%  x(1) = gaussian amplitude 1
%  x(2) = width 1 ( 1/(2*sigma^2) )
%  x(3) = centre freq of peak 1
%  x(4) = gaussian amplitude 2
%  x(5) = width 2 ( 1/(2*sigma^2) )
%  x(6) = centre freq of peak 2
%  x(7) = amplitude of linear baseline
%  x(8) = constant amplitude offset

% MM: Allowing peaks to vary individually seems to work better than keeping
% the distance fixed (i.e., including J in the function)

F = x(1)*exp(x(2)*(freq-x(3)).*(freq-x(3))) + ...
    x(4)*exp(x(5)*(freq-x(6)).*(freq-x(6))) + ...
    x(7)*(freq-x(3))+x(8);

%%%%%%%%%%%%%%%%%%%%%%%% DOUBLE GAUSS MODEL (MM: 150211) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = GABAGlxModel_area(x,freq)
%% Function for DoubleGaussModel Model -- with extra definitions


% A = 5/127.7; % for scalar coupling for Glu-C2 multiplet at 3.74 ppm (5 Hz seems to work best for in vivo data)
% F = x(1)*exp(x(2)*(freq-x(3)+A).*(freq-x(3)+A)) ...
%     +x(1)*exp(x(2)*(freq-x(3)-A).*(freq-x(3)-A)) ...
%     +x(4)*(freq-x(3))+x(5);

% Three Gaussians
%  x(1) = gaussian amplitude 1
%  x(2) = width 1 ( 1/(2*sigma^2) )
%  x(3) = centre freq of peak 1
%  x(4) = gaussian amplitude 2
%  x(5) = width 2 ( 1/(2*sigma^2) )
%  x(6) = centre freq of peak 2
%  x(7) = gaussian amplitude 3
%  x(8) = width 3 ( 1/(2*sigma^2) )
%  x(9) = centre freq of peak 3
%  x(10) = amplitude of linear baseline
%  x(11) = constant amplitude offset
%  x(12) = sine term
%  x(13) = cosine term



% MM: Allowing peaks to vary individually seems to work better than keeping
% the distance fixed (i.e., including J in the function)

F = x(1)*exp(x(2)*(freq-x(3)).*(freq-x(3))) + ...
    x(4)*exp(x(5)*(freq-x(6)).*(freq-x(6))) + ...
    x(7)*exp(x(8)*(freq-x(9)).*(freq-x(9))) + ...
    x(10)*(freq-x(3))+...
    x(11)*sin(pi*freq/1.31/4)+x(12)*cos(pi*freq/1.31/4);




%%%%%%%%%%%%%%% BASELINE %%%%%%%%%%%%%%%%%%%%%%%
function F = BaselineModel(x,freq)
%% Function for Baseline Model  

F = x(2)*(freq-x(1))+x(3);




%%%%%%%%%%%%%%% BASELINE %%%%%%%%%%%%%%%%%%%%%%%
function data_corr = klose_eddy_correction(data,H2O_fid)
%% Phase correction using Klose method

K=abs(H2O_fid);
Kphase=phase(H2O_fid);


K_data=abs(data);
Kphase_data=phase(data);

Kphase_corr=Kphase_data-Kphase;

data_corr=K_data.* exp(i*Kphase_corr);




%%%%%%%%%%%%%%%%%%% GABA INST UNITS CALC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MRS_struct] = MRSGABAinstunits(MRS_struct,reg,metab,ii)
%% Function for GABA Concentration estimation   

% function [MRS_struct] = MRSGABAinstunits(MRS_struct)
% Convert GABA and Water amplitudes to institutional units
% (pseudo-concentration in mmol per litre).
% March 10: use MRS_struct.
TR=MRS_struct.p.TR/1000;
TE=MRS_struct.p.TE/1000;
PureWaterConc = 55000; % mmol/litre
WaterVisibility = 0.65; % This is approx the value from Ernst, Kreis, Ross
T1_Water = 1.100; % average of WM and GM, estimated from Wansapura 1999
T2_Water = 0.095; % average of WM and GM, estimated from Wansapura 1999
N_H_Water=2;

switch metab
            case 'GABA'
                EditingEfficiency = 0.5;
                T1 = 0.80 ; % "empirically determined"...! Gives same values as RE's spreadsheet
                % ... and consistent with Cr-CH2 T1 of 0.8 (Traber, 2004)
                %Not yet putting in measured GABA T1, but it is in the pipeline - 1.35ish
                T2 = 0.13; % from occipital Cr-CH2, Traber 2004
                T2 = 0.088; % from JMRI paper 2012 Edden et al.
                N_H_metab=2;
                MM=0.45;  % MM correction: fraction of GABA in GABA+ peak. (In TrypDep, 30 subjects: 55% of GABA+ was MM)
                %This fraction is platform and implementation dependent, base on length and
                %shape of editing pulses and ifis Henry method.
                T1_factor = (1-exp(-TR./T1_Water)) ./ (1-exp(-TR./T1));
                T2_factor = exp(-TE./T2_Water) ./ exp(-TE./T2);
                
            case 'Glx'
                EditingEfficiency = 1.0;
                T1 = 1 ; % Not yet in the literature -- MGSaleh
                T2 = 1; % Not yet in the literature -- MGSaleh
                N_H_metab=2; % maybe -- MGSaleh                
                MM=1;  % Not affected by MM -- MGSaleh
                T1_factor = (1-exp(-TR./T1_Water)) ./ (1-exp(-TR./T1));
                T2_factor = exp(-TE./T2_Water) ./ exp(-TE./T2);
             
            case 'GSH'
                EditingEfficiency = 0.74;  %At 3T based on Quantification of Glutathione in the Human Brain by MR Spectroscopy at 3 Tesla: 
                %Comparison of PRESS and MEGA-PRESS
                %Faezeh Sanaei Nezhad etal. DOI 10.1002/mrm.26532, 2016 -- MGSaleh
                T1 = 0.40 ; % At 3T based on Doubly selective multiple quantum chemical shift imaging and 
                % T1 relaxation time measurement of glutathione (GSH) in the human brain in vivo
                % In-Young Choi et al. NMR Biomed. 2013; 26: 28?34 -- MGSaleh
                T2 = 0.12;  % At 3T based on the ISMRM abstract
                %T2 relaxation times of 18 brain metabolites determined in 83 healthy volunteers in vivo
                % Milan Scheidegger et al. Proc. Intl. Soc. Mag. Reson. Med. 22 (2014)-- MGSaleh
                N_H_metab=2; % Need to check in the literature -- MGSaleh
                MM= 1;
                T1_factor = (1-exp(-TR./T1_Water)) ./ (1-exp(-TR./T1));
                T2_factor = exp(-TE./T2_Water) ./ exp(-TE./T2);
                
            case 'Lac'
                EditingEfficiency = 1.0;
                T1 = 1 ; % Need to check in the literature -- MGSaleh
                T2 = 1; % Need to check in the literature -- MGSaleh
                N_H_metab=2; % Need to check in the literature -- MGSaleh
                MM=1;  % Not affected by MM -- MGSaleh
                T1_factor = (1-exp(-TR./T1_Water)) ./ (1-exp(-TR./T1));
                T2_factor = exp(-TE./T2_Water) ./ exp(-TE./T2);
end



Nspectra = length(MRS_struct.gabafile);
%Nwateravg=8;


if(strcmpi(MRS_struct.p.vendor,'Siemens'))
    MRS_struct.out.(reg).(metab).conciu(ii) = (MRS_struct.out.(reg).(metab).Area(ii)  ./  MRS_struct.out.(reg).Water.Area(ii))  ...
    * PureWaterConc*WaterVisibility*T1_factor*T2_factor*(N_H_Water./N_H_metab) ...
    * MM /2.0 ./ EditingEfficiency; %Factor of 2.0 is appropriate for averaged data, read in separately as on and off (Siemens).
else
    MRS_struct.out.(reg).(metab).conciu(ii) = (MRS_struct.out.(reg).(metab).Area(ii)  ./  MRS_struct.out.(reg).Water.Area(ii))  ...
    * PureWaterConc*WaterVisibility*T1_factor*T2_factor*(N_H_Water./N_H_metab) ...
    * MM ./ EditingEfficiency;
end
FAC=PureWaterConc*WaterVisibility*(N_H_Water./N_H_metab) ...
    * MM ./ EditingEfficiency*T1_factor*T2_factor;



%%%%%%%%%%%%%%%%%%% Glx INST UNITS CALC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [MRS_struct] = MRSGlxinstunits(MRS_struct,reg,metab,ii)
% %% Function for Glx Concentration estimation   
% 
% % function [MRS_struct] = MRSGlxinstunits(MRS_struct)
% % Convert Glx and Water amplitudes to institutional units
% % (pseudo-concentration in mmol per litre).
% %  Written by MGSaleh -- 29 June 2016. Still needs variable inputs, such as
% %  T1, etc.
% 
% PureWaterConc = 55000; % mmol/litre
% WaterVisibility = 0.65; % This is approx the value from Ernst, Kreis, Ross
% EditingEfficiency = 1; % Need to verify -- MGSaleh
% 
% T1_Glx = 1 ; % Not yet in the literature -- MGSaleh
% 
% T2_Glx = 1; % Not yet in the literature -- MGSaleh
% T2_Glx = 1; % Not yet in the literature -- MGSaleh
% 
% T1_Water = 1.100; % average of WM and GM, estimated from Wansapura 1999
% T2_Water = 0.095; % average of WM and GM, estimated from Wansapura 1999
% %
% TR=MRS_struct.p.TR/1000;
% TE=MRS_struct.p.TE/1000;
% 
% N_H_Water=2;
% Nspectra = length(MRS_struct.gabafile);
% %Nwateravg=8;
% 
% T1_factor = (1-exp(-TR./T1_Water)) ./ (1-exp(-TR./T1_Glx));
% T2_factor = exp(-TE./T2_Water) ./ exp(-TE./T2_Glx);
% 
% if(strcmpi(MRS_struct.p.vendor,'Siemens'))
%     MRS_struct.out.(reg).(metab).conciu(ii) = (MRS_struct.out.(reg).(metab).Area(ii)  ./  MRS_struct.out.(reg).Water.Area(ii))  ...
%     * PureWaterConc*WaterVisibility*T1_factor*T2_factor*(N_H_Water./N_H_Glx) ...
%     * MM /2.0 ./ EditingEfficiency; %Factor of 2.0 is appropriate for averaged data, read in separately as on and off (Siemens).
% else
%     MRS_struct.out.(reg).(metab).conciu(ii) = (MRS_struct.out.(reg).(metab).Area(ii)  ./  MRS_struct.out.(reg).Water.Area(ii))  ...
%     * PureWaterConc*WaterVisibility*T1_factor*T2_factor*(N_H_Water./N_H_Glx) ...
%     * MM ./ EditingEfficiency;
% end
% FAC=PureWaterConc*WaterVisibility*(N_H_Water./N_H_Glx) ...
%     * MM ./ EditingEfficiency*T1_factor*T2_factor;
% 

%%%%%%%%%%%%%%% INSET FIGURE %%%%%%%%%%%%%%%%%%%%%%%
function [h_main, h_inset]=inset(main_handle, inset_handle,inset_size)
%% Function for figure settings 


% The function plotting figure inside figure (main and inset) from 2 existing figures.
% inset_size is the fraction of inset-figure size, default value is 0.35
% The outputs are the axes-handles of both.
%
% An examle can found in the file: inset_example.m
%
% Moshe Lindner, August 2010 (C).

if nargin==2
    inset_size=0.35;
end

inset_size=inset_size*.5;
%figure
new_fig=gcf;
main_fig = findobj(main_handle,'Type','axes');
h_main = copyobj(main_fig,new_fig);
set(h_main,'Position',get(main_fig,'Position'))
inset_fig = findobj(inset_handle,'Type','axes');
h_inset = copyobj(inset_fig,new_fig);
ax=get(main_fig,'Position');
set(h_inset,'Position', [1.3*ax(1)+ax(3)-inset_size 1.001*ax(2)+ax(4)-inset_size inset_size*0.7 inset_size*0.9])







