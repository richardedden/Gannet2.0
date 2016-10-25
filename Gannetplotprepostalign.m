function Gannetplotprepostalign(MRS_struct, specno)
%function MRSplotprepostalign(MRS_struct, specno)
% Plots pre and post alignment spectra in MRSLoadPfiles
% 110214:  Scale spectra by the peak _height_ of water
%          Plot multiple spectra as a stack - baselines offset
%            by mean height of GABA

if MRS_struct.p.HERMES
        numspec = 4;
    
    
    
        SpectraToPlot = [eval(['MRS_struct.spec.', sprintf('%s',MRS_struct.p.target),'.diff']); ...
        eval(['MRS_struct.spec.', sprintf('%s',MRS_struct.p.target),'.diff_noalign']); ...
        eval(['MRS_struct.spec.', sprintf('%s',MRS_struct.p.target2),'.diff']); ...
        eval(['MRS_struct.spec.', sprintf('%s',MRS_struct.p.target2),'.diff_noalign']);];
    
        % Estimate baseline from between GABAGlx or Lac and GSH. The values might be changed depending on the future choice of metabolites -- MGSaleh
        
        if strcmp(MRS_struct.p.target2, 'Lac') % Defining different limits for diferent target -- MGSaleh 2016
            
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
    
        
        %To determine the output depending on the type of acquistion used -- MGSaleh 2016
        SpectraToPlot = [eval(['MRS_struct.spec.', sprintf('%s',MRS_struct.p.target),'.diff']); ...
        eval(['MRS_struct.spec.', sprintf('%s',MRS_struct.p.target),'.diff_noalign']);];
    
        % Estimate baseline from between Glx and GABA
        z=abs(MRS_struct.spec.freq-3.6);
        Glx_right=find(min(z)==z);
        z=abs(MRS_struct.spec.freq-3.3);
        GABA_left=find(min(z)==z);
        z=abs(MRS_struct.spec.freq-2.8);
        GABA_right=find(min(z)==z);
        specbaseline = (mean(real(SpectraToPlot(1,Glx_right:GABA_left)),2));
    
    
end


% SpectraToPlot = [MRS_struct.spec.diff(specno,:); MRS_struct.spec.diff_noalign(specno,:)];


% % Estimate baseline from between Glx and GABA
% z=abs(MRS_struct.spec.freq-3.6);
% Glx_right=find(min(z)==z);
% z=abs(MRS_struct.spec.freq-3.3);
% GABA_left=find(min(z)==z);
% z=abs(MRS_struct.spec.freq-2.8);
% GABA_right=find(min(z)==z);;
% specbaseline = (mean(real(SpectraToPlot(1,Glx_right:GABA_left)),2));





% Added by MGSaleh 2016
if MRS_struct.p.HERMES
    
    % averaged gaba height across all scans - to estimate stack spacing
    gabaheight = abs(max(SpectraToPlot(1,Glx_right:GABA_right),[],2));
    gabaheight = mean(gabaheight);
    plotstackoffset = [ 0 : (numspec-1) ]';
        
    if strcmp(MRS_struct.p.target2, 'Lac')  % Defining different limits for diferent target -- MGSaleh 2016
        plotstackoffset = plotstackoffset * 0.5 * gabaheight;
    else
        plotstackoffset = plotstackoffset * 1.0 * gabaheight;
    end
    plotstackoffset = plotstackoffset - specbaseline;

    
    plot(MRS_struct.spec.freq, real(SpectraToPlot((1),:)),'b',MRS_struct.spec.freq, real(SpectraToPlot((2),:)),'r');
    hold on
    shft=repmat(plotstackoffset, [1 length(SpectraToPlot(1,:))]);
    SpectraToPlot(3:4,:) = SpectraToPlot(3:4,:) + [max(shft,[],1); max(shft,[],1)] ;
    plot(MRS_struct.spec.freq, real(SpectraToPlot((3),:)),'b',MRS_struct.spec.freq, real(SpectraToPlot((4),:)),'r');
    hold off
    
    % yaxismax
    % yaxismin
    if strcmp(MRS_struct.p.target2, 'Lac')
        yaxismax = (numspec + 1.0) * 0.5 * gabaheight; % top spec + 0.5 * height of gaba %Changed slightly by MGSaleh to accomodate both GSH and GSH/Lac -- 2016
    else
        yaxismax = (numspec + 1.0) * 1.0 * gabaheight; % top spec + 1.0 * height of gaba %Changed slightly by MGSaleh to accomodate both GSH and GABAGlx/GSH -- 2016
    end
    yaxismin =  - 2.0* gabaheight; % extend 2* gaba heights below zero %Changed slightly by MGSaleh to accomodate both GSH and GABAGlx/Lac -- 2016
    
    if (yaxismax<yaxismin)
        dummy=yaxismin;
        yaxismin=yaxismax;
        yaxismax=dummy;
    end

    
else
    % averaged gaba height across all scans - to estimate stack spacing
    gabaheight = abs(max(SpectraToPlot(1,Glx_right:GABA_right),[],2));
    gabaheight = mean(gabaheight);
    plotstackoffset = [ 0 : (numspec-1) ]';
    plotstackoffset = plotstackoffset * gabaheight;
    plotstackoffset = plotstackoffset - specbaseline;
    
    SpectraToPlot = SpectraToPlot + ...
        repmat(plotstackoffset, [ 1  length(SpectraToPlot(1,:))]);
    
    plot(MRS_struct.spec.freq, real(SpectraToPlot));
    
    % yaxismax
    % yaxismin
    yaxismax = (numspec + 1.0) *gabaheight; % top spec + 2* height of gaba %Changed slightly by MGSaleh to accomodate both GSH and GABAGlx/Lac -- 2016
    yaxismin =  - 2.0* gabaheight; % extend 2* gaba heights below zero %Changed slightly by MGSaleh to accomodate both GSH and GABAGlx/Lac -- 2016
    if (yaxismax<yaxismin)
        dummy=yaxismin;
        yaxismin=yaxismax;
        yaxismax=dummy;
    end

end



legendtxt = {'post', 'pre'};
hl=legend(legendtxt);
set(hl,'EdgeColor',[1 1 1]);
set(gca,'XDir','reverse');
oldaxis = axis;



axis([0 5  yaxismin yaxismax])
