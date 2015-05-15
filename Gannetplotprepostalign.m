function Gannetplotprepostalign(MRS_struct, specno)
%function MRSplotprepostalign(MRS_struct, specno)
% Plots pre and post alignment spectra in MRSLoadPfiles
% 110214:  Scale spectra by the peak _height_ of water
%          Plot multiple spectra as a stack - baselines offset
%            by mean height of GABA

numspec = 2;
SpectraToPlot = [MRS_struct.spec.diff(specno,:); MRS_struct.spec.diff_noalign(specno,:)];


% Estimate baseline from between Glx and GABA
z=abs(MRS_struct.spec.freq-3.6);
Glx_right=find(min(z)==z);
z=abs(MRS_struct.spec.freq-3.3);
GABA_left=find(min(z)==z);
z=abs(MRS_struct.spec.freq-2.8);
GABA_right=find(min(z)==z);
specbaseline = mean(real(SpectraToPlot(1,Glx_right:GABA_left)),2);
% averaged gaba height across all scans - to estimate stack spacing
gabaheight = abs(max(SpectraToPlot(1,Glx_right:GABA_right),[],2));
gabaheight = mean(gabaheight);
plotstackoffset = [ 0 : (numspec-1) ]';
plotstackoffset = plotstackoffset * gabaheight;
plotstackoffset = plotstackoffset - specbaseline;
SpectraToPlot = SpectraToPlot + ...
    repmat(plotstackoffset, [ 1  length(MRS_struct.spec.diff(1,:))]);
plot(MRS_struct.spec.freq, real(SpectraToPlot));
legendtxt = {'post', 'pre'};
hl=legend(legendtxt);
set(hl,'EdgeColor',[1 1 1]);
set(gca,'XDir','reverse');
oldaxis = axis;
yaxismax = (numspec + 2) *gabaheight; % top spec + 2* height of gaba
yaxismin =  - 2* gabaheight; % extend 2* gaba heights below zero
if (yaxismax<yaxismin)
    dummy=yaxismin;
    yaxismin=yaxismax;
    yaxismax=dummy;
end
% yaxismax
% yaxismin
axis([0 5  yaxismin yaxismax])
