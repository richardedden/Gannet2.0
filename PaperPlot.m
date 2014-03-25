function PaperPlot(x,y,options)
%PaperPlot(MRS_struct.spec.freq,MRS_struct.spec.freq,'k')
%This plots a spectrum (or all spectra) from the input arguments.
%For paper output, save as .eps format (matlab pdf isn't good).

if nargin == 2
    plot(x,y)
else
    plot(x,y,options)
end   

set(gca,'XLim',[0.5 4.5]);
set(gca,'XDir','reverse');
set(gca,'YTick',[]);
set(gca,'Box','off');
