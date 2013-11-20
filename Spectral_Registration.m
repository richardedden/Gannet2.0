function [AllFramesFTrealign MRS_struct] = Spectral_Registration(MRS_struct)
    %Spectral Registration is a time-domain frequency-and-phse correction as
    %per Near et al. 2014 [under review].

    %First, take the complex data and turn it into a real matrix
    MRS_struct.flatdata(:,1,:)=real(MRS_struct.data);
    MRS_struct.flatdata(:,2,:)=imag(MRS_struct.data);
    %Correct to a point 10% into the file (seems better that the actual beginning)
    AlignRow=ceil(size(MRS_struct.flatdata,3)/10);
    MRS_struct.flattarget=squeeze(MRS_struct.flatdata(:,:,AlignRow)); 
    
    
    %Time domain Frequency and Phase Correction
    %Preliminary to fitting:
    parsGuess=[0 0]; %initial freq and phase guess 
    parsFit = zeros([size(MRS_struct.data,2) 2]);
    input.dwelltime=1/MRS_struct.sw;
    time=((0:1:(MRS_struct.npoints-1)).'/MRS_struct.sw);
    %Fitting to determine frequency and phase corrections.
    for corrloop=1:size(MRS_struct.data,2)
        target=MRS_struct.flattarget(:);
        transient=squeeze(MRS_struct.flatdata(:,:,corrloop));
        input.data=transient(:);
        parsFit(corrloop,:)=nlinfit(input,target,@FreqPhaseShiftNest,parsGuess);
    end
    size(MRS_struct.data)
    %Applyng frequency and phase corrections.
    for corrloop=1:size(MRS_struct.data,2)
        MRS_struct.data_align(:,corrloop)=MRS_struct.data(:,corrloop).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));
    end
    FullData = MRS_struct.data_align;
    FullData = FullData.* repmat( (exp(-(time)*MRS_struct.LB*pi)), [1 size(MRS_struct.data,2)]);
    AllFramesFTrealign=fftshift(fft(FullData,MRS_struct.ZeroFillTo,1),1);

    %In FD, move Cr to 3.02 and get phase 'right' as opposed to 'consistent'.
    ChoCrFitLimLow=2.6;
    ChoCrFitLimHigh=3.6;           
    %Still need ranges for Creatine align plot
    z=abs(MRS_struct.freq-ChoCrFitLimHigh);
    cclb=find(min(z)==z);
    z=abs(MRS_struct.freq-ChoCrFitLimLow);
    ccub=find(min(z)==z);
    freqrange=MRS_struct.freq(cclb:ccub);
    %Do some detective work to figure out the initial parameters
    ChoCrMeanSpec = mean(AllFramesFTrealign(cclb:ccub,:),2);
    Baseline_offset=real(ChoCrMeanSpec(1)+ChoCrMeanSpec(end))/2;
    Width_estimate=0.05;%ppm
    Area_estimate=(max(real(ChoCrMeanSpec))-min(real(ChoCrMeanSpec)))*Width_estimate*4;
    ChoCr_initx = [ Area_estimate Width_estimate 3.02 0 Baseline_offset 0 1].*[1 (2*MRS_struct.LarmorFreq) (MRS_struct.LarmorFreq) (180/pi) 1 1 1];
       
    ChoCrMeanSpecFit = FitChoCr(freqrange, ChoCrMeanSpec, ChoCr_initx,MRS_struct.LarmorFreq);
    MRS_struct.ChoCrMeanSpecFit = ChoCrMeanSpecFit./[1 (2*MRS_struct.LarmorFreq) (MRS_struct.LarmorFreq) (180/pi) 1 1 1];
%     figure(5)
%     plot(freqrange,ChoCrMeanSpec,freqrange,TwoLorentzModel(ChoCrMeanSpecFit./[1 (2*MRS_struct.LarmorFreq) (MRS_struct.LarmorFreq) (180/pi) 1 1 1],freqrange))
    
    AllFramesFTrealign=AllFramesFTrealign*exp(1i*pi/180*(ChoCrMeanSpecFit(4)));%phase
        ChoCrFreqShift = ChoCrMeanSpecFit(3);
        ChoCrFreqShift = ChoCrFreqShift - 3.02*MRS_struct.LarmorFreq;
        ChoCrFreqShift = ChoCrFreqShift ./ (MRS_struct.LarmorFreq*(MRS_struct.freq(2) - MRS_struct.freq(1) ));
        ChoCrFreqShift_points = round(ChoCrFreqShift)
    AllFramesFTrealign=circshift(AllFramesFTrealign, [-ChoCrFreqShift_points 0]);%freq
        %MRS_struct.ChoCrFreqShift_pts(ii) = ChoCrFreqShift_points;
    
    %Fit just the Cr in the aligned mean spectrum to get CrFWHMHz
    CrFitLimLow=2.6;
    CrFitLimHigh=3.11;           
    %Still need ranges for Creatine align plot
    z=abs(MRS_struct.freq-CrFitLimHigh);
    clb=find(min(z)==z);
    z=abs(MRS_struct.freq-CrFitLimLow);
    cub=find(min(z)==z);
    freqrange=MRS_struct.freq(clb:cub);
    Cr_initx = [ Area_estimate Width_estimate 3.02 0 Baseline_offset 0 ].*[1 (2*MRS_struct.LarmorFreq) (MRS_struct.LarmorFreq) (180/pi) 1 1 ];
    CrMeanSpec = mean(AllFramesFTrealign(clb:cub,:),2);
    CrMeanSpecFit = FitCr(freqrange, CrMeanSpec, Cr_initx);
    
    %Some Output
    MRS_struct.FreqStdevHz(MRS_struct.ii)=std(parsFit(:,1),1);
    MRS_struct.CrFWHMHz(MRS_struct.ii)=CrMeanSpecFit(2);


 
                
                
                


            
            
                
%             AllFramesFTrealign=reshape(AllFramesFTrealign,[size(AllFramesFTrealign,1) MRS_struct.Navg(ii) MRS_struct.nrows]);
% 
%             OddFramesFTrealign=AllFramesFTrealign(:,:,2:2:end);
%             EvenFramesFTrealign=AllFramesFTrealign(:,:,1:2:end);
%            
%             OddFramesFTrealign=reshape(OddFramesFTrealign,[size(OddFramesFTrealign,1) size(OddFramesFTrealign,2)*size(OddFramesFTrealign,3) ]);
%             EvenFramesFTrealign=reshape(EvenFramesFTrealign,[size(EvenFramesFTrealign,1) size(EvenFramesFTrealign,2)*size(EvenFramesFTrealign,3) ]);
%            
%             MRS_struct.Navg(ii) = MRS_struct.Navg(ii)*MRS_struct.nrows; % ADH fix to have the correct Navg for .data now    
              




end