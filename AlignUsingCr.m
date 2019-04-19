function AllFramesFTrealign=AlignUsingCr(AllFramesFTrealign,MRS_struct)
    %Using the freq information from MRS_struct, determine appropriate range to
    %fit Cr signal over
    Cr_ppm=3.02;
    %Determine limits
    z=abs(MRS_struct.spec.freq-3.12);
    lb=find(min(z)==z);
    z=abs(MRS_struct.spec.freq-2.72);
    ub=find(min(z)==z);
    %Set initial parameters by fitting the sum
    
    freqrange = MRS_struct.spec.freq(lb:ub);
    Initx = [ 30 0.05 3.02 0 0 0 ];
    SumSpec = sum(AllFramesFTrealign(lb:ub,:),2);
    SumSpecParams = FitPeaksByFrames2(freqrange, SumSpec, Initx); %FitPeaksByFrames2 now edited so that output params have same units as Initx
    %Update freq and phase of AllFramesFTrealign to reflect Sum fit
    AllFramesFTrealign=AllFramesFTrealign*exp(1i*SumSpecParams(4));
    %Shift Cr freq to 3.02
    freq_step_size=abs(MRS_struct.spec.freq(1)-MRS_struct.spec.freq(2));
    Shift_points=round((real(SumSpecParams(3))-Cr_ppm)/freq_step_size);
    AllFramesFTrealign=circshift(AllFramesFTrealign,[Shift_points 0]);
    %Set initial fitting parameters by fitting the sum
    Init = SumSpecParams./[ size(AllFramesFTrealign,2) 1 1 1 size(AllFramesFTrealign,2) size(AllFramesFTrealign,2) ];
    Init(3)=Cr_ppm;
    Init(4)=0;
    %Select Data2BeFit
    %subplot(2,2,1)
    %plot(MRS_struct.spec.freq,real(AllFramesFTrealign))
    %set(gca,'XLim',[2.5 3.5]);
    %set(gca,'XDir','reverse')
    Data2bFit = AllFramesFTrealign(lb:ub,:);
    FrameParams = FitPeaksByFrames2(freqrange, Data2bFit, Init); %FitPeaksByFrames2 now edited so that output params have same units as Initx
    %Update freq and phase of AllFramesFTrealign to reflect Frame-by-Frame fit
    FrameShift_points=round((real(FrameParams(:,3))-Cr_ppm)/freq_step_size);
    for jj=1:size(FrameShift_points,1)
        AllFramesFTrealign(:,jj)=circshift(AllFramesFTrealign(:,jj),[-FrameShift_points(jj) 0]);
    end
    %subplot(2,2,2)
    %plot(MRS_struct.spec.freq,real(AllFramesFTrealign),'.-')
    %set(gca,'XLim',[2.5 3.5]);
    %set(gca,'XDir','reverse')
    %mean(FrameShift_points)
    %mean(FrameParams(:,3))
    %size(AllFramesFTrealign)
    %size(repmat(exp(1i*SumSpecParams(:,4)).',[length(MRS_struct.spec.freq) 1]))
    AllFramesFTrealign=AllFramesFTrealign.*repmat(exp(-1i*FrameParams(:,4)).',[length(MRS_struct.spec.freq) 1]);
    AllFramesFTrealign=AllFramesFTrealign+repmat(FrameParams(:,5).',[length(MRS_struct.spec.freq) 1]);
    %subplot(2,2,3)
    %plot(MRS_struct.spec.freq,real(AllFramesFTrealign))
    %set(gca,'XLim',[2.5 3.5]);
    %set(gca,'XDir','reverse');
    %plot the models.
    Models=zeros(size(Data2bFit));
    for ii=1:size(Models,2)
        Models(:,ii)=LorentzModel(FrameParams(ii,:),freqrange);
    end
    subplot(3,1,1)
    plot(freqrange,Models);
    set(gca,'XLim',[2.5 3.5]);
    set(gca,'XDir','reverse');
    for ii=1:size(Models,2)
    %Models_undone(:,ii)=((Models(:,ii)-FrameParams(ii,5)-FrameParams(ii,6)*(freqrange.'-FrameParams(ii,3))))*exp(-1i*FrameParams(ii,4));
    Models_undone(:,ii)=(Models(:,ii)*exp(1i*FrameParams(ii,4)));
    %Models_undone(:,ii)=circshift(Models_undone(:,ii),[FrameShift_points(jj) 0]);
    end
    subplot(3,1,2)
    plot(freqrange,Models_undone);
    set(gca,'XLim',[2.5 3.5]);
    set(gca,'XDir','reverse');
   subplot(3,1,3)
   plot(1:length(freqrange),std(Models.'),1:length(freqrange),std(Models_undone.'))
   plot(freqrange,sum(Data2bFit,2),freqrange,sum(Models_undone,2))
    save testing.mat Data2bFit FrameParams freqrange
    %plot(FrameParams(:,4)/pi*180)
    
end