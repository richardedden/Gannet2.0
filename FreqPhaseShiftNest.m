function output=FreqPhaseShiftNest(pars,input);
        %Heavily copied from jn_freqPhaseShiftNest of Jamie Near, McGill.
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        dwelltime=input.dwelltime;
        t=[0:dwelltime:(length(input.data)/2-1)*dwelltime];
        input.data=reshape(input.data,[length(input.data)/2 2]);
        fid=input.data(:,1)+1i*input.data(:,2);
        
        y=fid.*exp(1i*pi*(t'*f*2+p/180));
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        output=[real(y) imag(y)];
        output=output(:);
    end

