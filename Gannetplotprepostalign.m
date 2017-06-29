function Gannetplotprepostalign(MRS_struct, reg, specno)
% Plots pre and post alignment spectra in MRSLoadPfiles
% 110214:  Scale spectra by the peak _height_ of water
%          Plot multiple spectra as a stack - baselines offset
%            by mean height of GABA
% Updates by MGSaleh 2016

for kk = 1:length(reg)
    
    if MRS_struct.p.HERMES
        
        numspec = 4;        
        SpectraToPlot = [MRS_struct.spec.(reg{kk}).(sprintf('%s',MRS_struct.p.target)).diff(specno,:); ...
            MRS_struct.spec.(reg{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(specno,:); ...
            MRS_struct.spec.(reg{kk}).(sprintf('%s',MRS_struct.p.target2)).diff(specno,:); ...
            MRS_struct.spec.(reg{kk}).(sprintf('%s',MRS_struct.p.target2)).diff_noalign(specno,:)];
        
        % Estimate baseline from between GABAGlx or Lac and GSH. The values might be changed depending on the future choice of metabolites
        if strcmp(MRS_struct.p.target2, 'Lac')         
            z=abs(MRS_struct.spec.freq-1.5);
            Glx_right=find(min(z)==z);
            z=abs(MRS_struct.spec.freq-1.0);
            GABA_left=find(min(z)==z);
            z=abs(MRS_struct.spec.freq-0.5);
            GABA_right=find(min(z)==z);            
        else
            z=abs(MRS_struct.spec.freq-3.1);
            Glx_right=find(min(z)==z);
            z=abs(MRS_struct.spec.freq-2.9);
            GABA_left=find(min(z)==z);
            z=abs(MRS_struct.spec.freq-2.8);
            GABA_right=find(min(z)==z);
        end        
        
        specbaseline = (mean(real(SpectraToPlot(1,Glx_right:GABA_left)),2));
        
    else
        
        numspec = 2;
        SpectraToPlot = [MRS_struct.spec.(reg{kk}).(sprintf('%s',MRS_struct.p.target)).diff(specno,:); ...
            MRS_struct.spec.(reg{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(specno,:)];
        
        % Estimate baseline from between Glx and GABA
        z=abs(MRS_struct.spec.freq-3.6);
        Glx_right=find(min(z)==z);
        z=abs(MRS_struct.spec.freq-3.3);
        GABA_left=find(min(z)==z);
        z=abs(MRS_struct.spec.freq-2.8);
        GABA_right=find(min(z)==z);
        specbaseline = (mean(real(SpectraToPlot(1,Glx_right:GABA_left)),2));
        
    end
    
    if MRS_struct.p.HERMES
        
        % averaged gaba height across all scans - to estimate stack spacing
        gabaheight = abs(max(SpectraToPlot(1,Glx_right:GABA_right),[],2));
        gabaheight = mean(gabaheight);
        plotstackoffset = (0:(numspec-1))';
        
        if strcmp(MRS_struct.p.target2, 'Lac')
            plotstackoffset = plotstackoffset * 0.5 * gabaheight;
        else
            plotstackoffset = plotstackoffset * 1.75 * gabaheight;
        end
        plotstackoffset = plotstackoffset - specbaseline;
        
        aa=1.2;
        plot(MRS_struct.spec.freq, aa*real(SpectraToPlot((2),:)), 'b', MRS_struct.spec.freq, aa*real(SpectraToPlot((1),:)), 'r');
        hold on;
        shift=repmat(plotstackoffset, [1 length(SpectraToPlot(1,:))]);
        SpectraToPlot(3:4,:) = SpectraToPlot(3:4,:) + [max(shift,[],1); max(shift,[],1)] ;
        plot(MRS_struct.spec.freq, aa*real(SpectraToPlot((4),:)),'b', MRS_struct.spec.freq, aa*real(SpectraToPlot((3),:)) ,'r');
        hold off;
        
        if strcmp(MRS_struct.p.target2, 'Lac')
            yaxismax = (numspec + 1.0) * 0.5 * gabaheight;
        else
            yaxismax = (numspec + 1.0) * 1.75 * gabaheight;
            
        end
        yaxismin = -2*gabaheight;        
        if yaxismax < yaxismin
            dummy=yaxismin;
            yaxismin=yaxismax;
            yaxismax=dummy;
        end        
        
    else
        % averaged gaba height across all scans - to estimate stack spacing
        gabaheight = abs(max(SpectraToPlot([1 2],Glx_right:GABA_right),[],2));
        gabaheight = max(gabaheight);
        plotstackoffset = (0:(numspec-1))';
        plotstackoffset = plotstackoffset * gabaheight;
        plotstackoffset = plotstackoffset - specbaseline;
        
        SpectraToPlot = SpectraToPlot + repmat(plotstackoffset, [1  length(SpectraToPlot(1,:))]);        
        plot(MRS_struct.spec.freq, real(SpectraToPlot(2,:)), 'b', MRS_struct.spec.freq, real(SpectraToPlot(1,:)), 'r');        
        
        yaxismax = 1.5*abs(max(max(real(SpectraToPlot([1 2],Glx_right:GABA_right)),[],2)));
        yaxismin = -10.0*abs(min(min(real(SpectraToPlot([1 2],Glx_right:GABA_right)),[],2)));
        if (yaxismax<yaxismin)
            dummy=yaxismin;
            yaxismin=yaxismax;
            yaxismax=dummy;
        end
        
    end
    
    legendtxt = {'pre', 'post'};
    hl = legend(legendtxt);
    set(hl,'EdgeColor',[1 1 1]);
    set(gca,'XDir','reverse');
    axis([0 5 yaxismin yaxismax]);
    
end
