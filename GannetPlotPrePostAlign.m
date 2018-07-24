function GannetPlotPrePostAlign(MRS_struct, vox, ii, kk)
% Plot pre-/post-alignment spectra
% 110214: Plot multiple spectra as a stack
% Updates by MGSaleh 2016, MM 2017-2018

if MRS_struct.p.HERMES
    
    SpectraToPlot = [MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:); ...
        MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:); ...
        MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff(ii,:); ...
        MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff_noalign(ii,:)];
    
    % Shift baselines to zero (MM: 180108)
    baserange = MRS_struct.spec.freq <= 0 & MRS_struct.spec.freq >= -0.5;
    switch MRS_struct.p.target2
        case {'GABA','GABAGlx','Glx'}
            if MRS_struct.p.phantom
                peakrange = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 1.0;
            else
                peakrange = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.79;
            end
        case 'GSH'
            peakrange = MRS_struct.spec.freq <= 3.5 & MRS_struct.spec.freq >= 1.5;
        case 'Lac'
            peakrange = MRS_struct.spec.freq <= 1.5 & MRS_struct.spec.freq >= 0.5;
    end
    
    specbaseline = mean(real(SpectraToPlot(:,baserange)),2);
    SpectraToPlot = SpectraToPlot - repmat(specbaseline, [1 length(SpectraToPlot)]);
    
    % Stack spectra (MM: 180724)
    if MRS_struct.p.phantom
        signalrange = max(max(real(SpectraToPlot(1:2,peakrange)))) - min(min(real(SpectraToPlot(1:2,peakrange))));
        SpectraToPlot(3:4,:) = SpectraToPlot(3:4,:) + signalrange;
    else
        peakheight = abs(min(real(SpectraToPlot(3,peakrange))));
        SpectraToPlot(3:4,:) = SpectraToPlot(3:4,:) + 1.05*peakheight;
    end
    hold on;
    plot(MRS_struct.spec.freq, real(SpectraToPlot(2,:)), 'Color', 'r');
    plot(MRS_struct.spec.freq, real(SpectraToPlot(1,:)), 'Color', 'b');
    plot(MRS_struct.spec.freq, real(SpectraToPlot(4,:)), 'Color', 'r');
    plot(MRS_struct.spec.freq, real(SpectraToPlot(3,:)), 'Color', 'b');
    hold off;
    
    yaxismax = abs(max(real(SpectraToPlot(3,peakrange))));
    yaxismax = yaxismax + yaxismax/4;
    yaxismin = abs(max(real(SpectraToPlot(1,peakrange))));
    yaxismin = -yaxismin/2;
    
else
    
    SpectraToPlot = [MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:); ...
        MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:)];
    
    % Shift baselines to zero (MM: 180108)
    baserange = MRS_struct.spec.freq <= 0 & MRS_struct.spec.freq >= -0.5;
    % Some bandwidth-limited acquisitions may not record anything below
    % 0 ppm, in this case get the baseline from the other side of
    % water. (GO: 180213)
    if sum(baserange) == 0
        baserange = MRS_struct.spec.freq >= 7 & MRS_struct.spec.freq <= 8;
    end
    switch MRS_struct.p.target
        case {'GABA','GABAGlx','Glx'}
            peakrange = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.79;
        case 'GSH'
            peakrange = MRS_struct.spec.freq <= 3.5 & MRS_struct.spec.freq >= 1.5;
        case 'Lac'
            peakrange = MRS_struct.spec.freq <= 4.5 & MRS_struct.spec.freq >= 0.5;
    end
    
    specbaseline = mean(real(SpectraToPlot(:,baserange)),2);
    SpectraToPlot = SpectraToPlot - repmat(specbaseline, [1 length(SpectraToPlot)]);
    
    % Stack spectra (MM: 180108)
    peakheight = max(abs(real(SpectraToPlot(1,peakrange))));
    SpectraToPlot(2,:) = SpectraToPlot(2,:) + peakheight;
    hold on;
    plot(MRS_struct.spec.freq, real(SpectraToPlot(2,:)), 'Color', 'r');
    plot(MRS_struct.spec.freq, real(SpectraToPlot(1,:)), 'Color', 'b');
    hold off;
    
    yaxismax = abs(max(real(SpectraToPlot(2,peakrange))));
    yaxismax = yaxismax + yaxismax/10;
    yaxismin = -peakheight/2;
    
end

box on;
legendtxt = {'pre','post'};
hl = legend(legendtxt);
set(hl,'EdgeColor',[1 1 1]);
set(gca,'XDir','reverse');
axis([0 5 yaxismin yaxismax]);



