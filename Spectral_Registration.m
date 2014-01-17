function [AllFramesFTrealign MRS_struct] = Spectral_Registration(MRS_struct, OnWhat, Dual)
    %Spectral Registration is a time-domain frequency-and-phse correction as
    %per Near et al. 2014 [under review].
    % OnWhat = 0 for spectro data OnWhat=1 for water data
    
 %Dual-channel option only applies registration separately to ONs and OFFs
    SpecRegLoop=0;
    if(nargin==3)
       if(Dual==1)
          %We want to run this code twice, once for ONs, once for OFFs.
          SpecRegLoop=1;
       end
    end
    while(SpecRegLoop>(-1))
    
    if(OnWhat) %Read water data
        %First, take the complex data and turn it into a real matrix
        flatdata(:,1,:)=real(MRS_struct.fids.data_water);
        flatdata(:,2,:)=imag(MRS_struct.fids.data_water);
        
        
    else % read spectro data
        
        if(nargin==3)
            if(Dual==1)
                %This code runs twice, first for ONs, second for OFFs.
                    SpecRegLoop;
                    size(real(MRS_struct.fids.data(:,(MRS_struct.fids.ON_OFF==SpecRegLoop))));
                    
                    flatdata(:,1,:)=real(MRS_struct.fids.data(:,(MRS_struct.fids.ON_OFF==SpecRegLoop)));
                    flatdata(:,2,:)=imag(MRS_struct.fids.data(:,(MRS_struct.fids.ON_OFF==SpecRegLoop)));
            end
        else
            %First, take the complex data and turn it into a real matrix
            flatdata(:,1,:)=real(MRS_struct.fids.data);
            flatdata(:,2,:)=imag(MRS_struct.fids.data);
        end  
    end
    
    %Correct to a point 10% into the file (seems better that the actual beginning)
    AlignRow=ceil(size(flatdata,3)/10);
    MRS_struct.fids.flattarget=squeeze(flatdata(:,:,AlignRow));
    
    
    %Time domain Frequency and Phase Correction
    %Preliminary to fitting:
    parsGuess=[0 0]; %initial freq and phase guess
    parsFit = zeros([size(flatdata,3) 2]);
    input.dwelltime=1/MRS_struct.p.sw;
    time=((0:1:(MRS_struct.p.npoints-1)).'/MRS_struct.p.sw);
    %Fitting to determine frequency and phase corrections.
    for corrloop=1:size(flatdata,3)
        target=MRS_struct.fids.flattarget(:);
        transient=squeeze(flatdata(:,:,corrloop));
        input.data=transient(:);
        parsFit(corrloop,:)=nlinfit(input,target,@FreqPhaseShiftNest,parsGuess);
    end   
    if(OnWhat)
        %Applyng frequency and phase corrections.
        for corrloop=1:size(flatdata,3)
            MRS_struct.fids.data_water(:,corrloop)=MRS_struct.fids.data_water(:,corrloop).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));
        end
        FullData = MRS_struct.fids.data_water;
        FullData = FullData.* repmat( (exp(-(time)*MRS_struct.p.LB*pi)), [1 size(MRS_struct.fids.data_water,2)]);
        AllFramesFTrealign=fftshift(fft(FullData,MRS_struct.p.ZeroFillTo,1),1);
        
    else
        
        %Applying frequency and phase corrections.
        for corrloop=1:size(flatdata,3)
            if(nargin==3)
                if(Dual==1)
                %Need to get the slot right to put data back into
                averages_per_dynamic=find(MRS_struct.fids.ON_OFF~=(MRS_struct.fids.ON_OFF(1)),1)-1;
                dyn=floor((corrloop-1)/averages_per_dynamic); %number of cycles in
                ind=mod((corrloop-1),averages_per_dynamic)+1; %number in current cycle
                
                if(SpecRegLoop==1)
                    if(MRS_struct.fids.ON_OFF(1)==1)
                        corrloop_d = dyn*averages_per_dynamic*2+ind;
                    else
                        corrloop_d = dyn*averages_per_dynamic*2+averages_per_dynamic+ind;
                    end
                else
                    if(MRS_struct.fids.ON_OFF(1)==1)
                        corrloop_d = dyn*averages_per_dynamic*2+averages_per_dynamic+ind;
                    else
                        corrloop_d = dyn*averages_per_dynamic*2+ind;
                    end
                end
                
                
                MRS_struct.fids.data_align(:,corrloop_d)=MRS_struct.fids.data(:,corrloop_d).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));       end
                CorrPars(corrloop_d,:)=parsFit(corrloop,:);
            else
                corrloop_d=corrloop;
                MRS_struct.fids.data_align(:,corrloop_d)=MRS_struct.fids.data(:,corrloop_d).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));
                CorrPars(corrloop_d,:)=parsFit(corrloop,:);
            end
            
        end
        
        if(SpecRegLoop==0)
            FullData = MRS_struct.fids.data_align;
            FullData = FullData.* repmat( (exp(-(time)*MRS_struct.p.LB*pi)), [1 size(MRS_struct.fids.data,2)]);
            AllFramesFTrealign=fftshift(fft(FullData,MRS_struct.p.ZeroFillTo,1),1);

            
            %In FD, move Cr to 3.02 and get phase 'right' as opposed to 'consistent'.
            ChoCrFitLimLow=2.6;
            ChoCrFitLimHigh=3.6;
            %Still need ranges for Creatine align plot
            z=abs(MRS_struct.spec.freq-ChoCrFitLimHigh);
            cclb=find(min(z)==z);
            z=abs(MRS_struct.spec.freq-ChoCrFitLimLow);
            ccub=find(min(z)==z);
            freqrange=MRS_struct.spec.freq(cclb:ccub);
            %Do some detective work to figure out the initial parameters
            ChoCrMeanSpec = mean(AllFramesFTrealign(cclb:ccub,:),2);
            Baseline_offset=real(ChoCrMeanSpec(1)+ChoCrMeanSpec(end))/2;
            Width_estimate=0.05;%ppm
            Area_estimate=(max(real(ChoCrMeanSpec))-min(real(ChoCrMeanSpec)))*Width_estimate*4;
            ChoCr_initx = [ Area_estimate Width_estimate 3.02 0 Baseline_offset 0 1].*[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1 1];

            
            if(nargin==3)
                if(Dual==1)
                    %This bit is silly - we don't want to do OFF-to-ON based on the CR signal      
                    ChoCrMeanSpecON = mean(AllFramesFTrealign(cclb:ccub,(MRS_struct.fids.ON_OFF==1)),2);
                    ChoCrMeanSpecOFF = mean(AllFramesFTrealign(cclb:ccub,(MRS_struct.fids.ON_OFF==0)),2);  
                    ChoCrMeanSpecFitON = FitChoCr(freqrange, ChoCrMeanSpecON, ChoCr_initx,MRS_struct.p.LarmorFreq);
                    ChoCrMeanSpecFitOFF = FitChoCr(freqrange, ChoCrMeanSpecOFF, ChoCr_initx,MRS_struct.p.LarmorFreq);
                    AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1))=AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1))*exp(1i*pi/180*(ChoCrMeanSpecFitON(4)));%phase
                    AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0))=AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0))*exp(1i*pi/180*(ChoCrMeanSpecFitOFF(4)));%phase
                    
                    ChoCrFreqShiftON = ChoCrMeanSpecFitON(3);
                    ChoCrFreqShiftON = ChoCrFreqShiftON - 3.02*MRS_struct.p.LarmorFreq;
                    ChoCrFreqShiftON = ChoCrFreqShiftON ./ (MRS_struct.p.LarmorFreq*(MRS_struct.spec.freq(2) - MRS_struct.spec.freq(1) ));
                    ChoCrFreqShift_pointsON = round(ChoCrFreqShiftON);
                    AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1))=circshift(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1)), [-ChoCrFreqShift_pointsON 0]);%freq
                    ChoCrFreqShiftOFF = ChoCrMeanSpecFitOFF(3);
                    ChoCrFreqShiftOFF = ChoCrFreqShiftOFF - 3.02*MRS_struct.p.LarmorFreq;
                    ChoCrFreqShiftOFF = ChoCrFreqShiftOFF ./ (MRS_struct.p.LarmorFreq*(MRS_struct.spec.freq(2) - MRS_struct.spec.freq(1) ));
                    ChoCrFreqShift_pointsOFF = round(ChoCrFreqShiftOFF);
                    AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0))=circshift(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0)), [-ChoCrFreqShift_pointsOFF 0]);%freq
            
                end
                
            MRS_struct.out.FreqStdevHz(MRS_struct.ii)=std(parsFit(:,1),1);
            MRS_struct.out.CrFWHMHz(MRS_struct.ii)=mean([ChoCrMeanSpecFitON(2) ChoCrMeanSpecFitOFF(2)]);    
            else
            ChoCrMeanSpecFit = FitChoCr(freqrange, ChoCrMeanSpec, ChoCr_initx,MRS_struct.p.LarmorFreq);
            MRS_struct.out.ChoCrMeanSpecFit = ChoCrMeanSpecFit./[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1 1];
            AllFramesFTrealign=AllFramesFTrealign*exp(1i*pi/180*(ChoCrMeanSpecFit(4)));%phase
            ChoCrFreqShift = ChoCrMeanSpecFit(3);
            ChoCrFreqShift = ChoCrFreqShift - 3.02*MRS_struct.p.LarmorFreq;
            ChoCrFreqShift = ChoCrFreqShift ./ (MRS_struct.p.LarmorFreq*(MRS_struct.spec.freq(2) - MRS_struct.spec.freq(1) ));
            ChoCrFreqShift_points = round(ChoCrFreqShift);
            AllFramesFTrealign=circshift(AllFramesFTrealign, [-ChoCrFreqShift_points 0]);%freq
            end
            %     figure(5)
            %     plot(freqrange,ChoCrMeanSpec,freqrange,TwoLorentzModel(ChoCrMeanSpecFit./[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1 1],freqrange))

 
            %Fit just the Cr in the aligned mean spectrum to get CrFWHMHz
            CrFitLimLow=2.6;
            CrFitLimHigh=3.11;
            %Still need ranges for Creatine align plot
            z=abs(MRS_struct.spec.freq-CrFitLimHigh);
            clb=find(min(z)==z);
            z=abs(MRS_struct.spec.freq-CrFitLimLow);
            cub=find(min(z)==z);
            freqrange=MRS_struct.spec.freq(clb:cub);
            Cr_initx = [ Area_estimate Width_estimate 3.02 0 Baseline_offset 0 ].*[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1 ];
            CrMeanSpec = mean(AllFramesFTrealign(clb:cub,:),2);
            CrMeanSpecFit = FitCr(freqrange, CrMeanSpec, Cr_initx);
            
            %Some Output
            MRS_struct.out.FreqStdevHz(MRS_struct.ii)=std(parsFit(:,1),1);
            MRS_struct.out.CrFWHMHz(MRS_struct.ii)=CrMeanSpecFit(2);
            %Decide which rows to reject based on 3-sigma
            
            
            
            % Reject any point where the fit params - phase
            %  or freq are > 3stdev away from the mean
            MeanFrameParams = mean(CorrPars, 1);
            UpperLim = repmat(MeanFrameParams + 3*std(CorrPars,1), [size(AllFramesFTrealign,2) 1]);
            LowerLim = repmat(MeanFrameParams - 3*std(CorrPars,1), [size(AllFramesFTrealign,2) 1]);
            %but don't reject on linear, const baseline fit vals
            rejectframe = gt(CorrPars, UpperLim);
            rejectframe = rejectframe + lt(CorrPars, LowerLim);
            MRS_struct.out.reject(:,MRS_struct.ii) = max(rejectframe,[],2);
            %Balance up rejects
            IsRejectOnOrOff=MRS_struct.fids.ON_OFF(MRS_struct.out.reject(:,MRS_struct.ii)==1);
            stepsize=find(MRS_struct.fids.ON_OFF~=(MRS_struct.fids.ON_OFF(1)),1)-1;
            for kk=1:size(MRS_struct.out.reject(:,MRS_struct.ii),1)
                %first find whether reject is ON or OFF
                if(MRS_struct.out.reject(kk,MRS_struct.ii)==1)
                    IsRejectOnOrOff=MRS_struct.fids.ON_OFF(kk);
                    if IsRejectOnOrOff==MRS_struct.fids.ON_OFF(1)
                         MRS_struct.out.reject(kk+stepsize,MRS_struct.ii)=1;
                    else
                        MRS_struct.out.reject(kk-stepsize,MRS_struct.ii)=1;
                    end
                end
            end
            
            
            
            
            
            
        end
    end
        SpecRegLoop=SpecRegLoop-1;
    end
    
  

end