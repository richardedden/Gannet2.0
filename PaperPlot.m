function PaperPlot(x,y,options)
%PaperPlot(MRS_struct) or PaperPlot(MRS_struct.spec.freq,MRS_struct.spec.diff,'k')

%This plots a spectrum (or all spectra) from the input arguments.
%For paper output, save as .eps format (matlab pdf isn't good).

switch nargin
    case 1
        z=abs(x.spec.freq-3.55);
        lowerbound=find(min(z)==z);
        z=abs(x.spec.freq-2.79);%2.75
        upperbound=find(min(z)==z);
        freqbounds=lowerbound:upperbound;
        freq=x.spec.freq(1,freqbounds);
        plot(x.spec.freq(1,:),x.spec.diff(11,:),'k',freq,GaussModel(x.out.GABAModelFit(11,:),freq),'r');
    case 2
    plot(x,y)
    case 3
    plot(x,y,options)
end   

set(gca,'XLim',[0.5 4.5]);
set(gca,'XDir','reverse');
set(gca,'YTick',[]);
set(gca,'Box','off');
