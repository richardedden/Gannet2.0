function MRS_struct = GannetQuantify(MRS_struct)

MRS_struct.version.quantify = '170831';

cWM = 1; % concentration of GABA in pure WM
cGM = 2; % concentration of GABA in pure GM
% NB: cWM/cGM is term used ratio not absolute value is important piece
alpha = cWM/cGM;

% Constants
% From Wansapura 1999  JMRI; 9:531 (1999)
%        T1          T2
% WM   832 +/- 10  79.2 +/- 0.6
% GM  1331 +/- 13  110 +/- 2
%
% From Lu, JMRI; 2005; 22: 13
% CSF T1 = 3817 +/- 424msec - but state may underestimated and that 4300ms
% is likely more accurate - but the reference is to an ISMRM 2001 abstract
% MacKay (last author) 2006 ISMRM abstract has T1 CSF = 3300 ms
% CSF T2 = 503.0 +/- 64.3 Piechnik MRM 2009; 61: 579
% However, other values from Stanisz et al:
% CPMG for T2, IR for T1
% T2GM = 99 +/ 7, lit: 71+/- 10 (27)
% T1GM = 1820 +/- 114, lit 1470 +/- 50 (29)
% T2WM = 69 +/-3 lit 56 +/- 4 (27)
% T1WM = 1084 +/- 45 lit 1110 +/- 45 (29)

T1w_WM    = 0.832;
T2w_WM    = 0.0792;
T1w_GM    = 1.331;
T2w_GM    = 0.110;
T1w_CSF   = 3.817;
T2w_CSF   = 0.503;
N_H_Water = 2;

% Determine concentration of water in GM, WM and CSF
% Gasparovic et al, MRM 2006; 55:1219 uses relative densities, ref to Ernst
% fGM = 0.78
% fWM = 0.65
% fCSH = 0.97
% such that
% concw_GM = 0.78 * 55.51 mol/kg = 43.30
% concw_WM = 0.65 * 55.51 mol/kg = 36.08
% concw_CSF = 0.97 * 55.51 mol/kg = 53.84

concw_GM  = 43.30*1e3;
concw_WM  = 36.08*1e3;
concw_CSF = 53.84*1e3;


if MRS_struct.p.PRIAM
    vox = {MRS_struct.p.Vox};
else
    vox = {MRS_struct.p.Vox{1}};
end

if MRS_struct.p.HERMES
    target = {MRS_struct.p.target, MRS_struct.p.target2};
else
    target = {MRS_struct.p.target};
end

tmp = strcmp(target,'GABAGlx');
if any(tmp)
    if MRS_struct.p.HERMES
        target = {'GABA','Glx',target{~tmp}};
    else
        target = {'GABA','Glx'};
    end
end

for ii = 1:length(MRS_struct.metabfile)
    
    TR = MRS_struct.p.TR(ii)/1000;
    TE = MRS_struct.p.TE(ii)/1000;
    
    for kk = 1:length(vox)
        
        meanGMfra = mean(MRS_struct.out.(vox{kk}).tissue.GMfra); % fraction of voxel that is GM for voxel fraction normalization - move to preinitialize
        meanWMfra = mean(MRS_struct.out.(vox{kk}).tissue.WMfra);
        
        fracGM  = MRS_struct.out.(vox{kk}).tissue.GMfra(ii);
        fracWM  = MRS_struct.out.(vox{kk}).tissue.WMfra(ii);
        fracCSF = MRS_struct.out.(vox{kk}).tissue.CSFfra(ii);
        
        CorrFactor = (meanGMfra + alpha*meanWMfra) / ((fracGM + alpha*fracWM) * (meanGMfra + meanWMfra));
        
        for trg = 1:length(target)
            
            switch target{trg}
                case 'GABA'
                    EditingEfficiency = 0.5; % For TE = 68 ms
                    T1_Metab = 1.31;  % Puts et al. 2013 (JMRI)
                    T2_Metab = 0.088; % Edden et al. 2012 (JMRI)
                    N_H_Metab = 2;
                    MM = 0.45; % MM correction: fraction of GABA in GABA+ peak. (In TrypDep, 30 subjects: 55% of GABA+ was MM)
                    % This fraction is platform and implementation dependent, based on length and
                    % shape of editing pulses and ifis Henry method
                    
                case 'Glx'
                    EditingEfficiency = 0.4; % determined by FID-A simulations (for TE = 68 ms)
                    T1_Metab = 1.23; % Posse et al. 2007 (MRM)
                    T2_Metab = 0.18; % Ganji et al. 2012 (NMR Biomed)
                    N_H_Metab = 1;
                    MM = 1;
                    
                case 'GSH'
                    EditingEfficiency = 0.74;  % At 3T based on Quantification of Glutathione in the Human Brain by MR Spectroscopy at 3 Tesla:
                    % Comparison of PRESS and MEGA-PRESS
                    % Faezeh Sanaei Nezhad etal. DOI 10.1002/mrm.26532, 2016 -- MGSaleh
                    T1_Metab = 0.40 ; % At 3T based on Doubly selective multiple quantum chemical shift imaging and
                    % T1 relaxation time measurement of glutathione (GSH) in the human brain in vivo
                    % In-Young Choi et al. NMR Biomed. 2013; 26: 28?34 -- MGSaleh
                    T2_Metab = 0.12; % At 3T based on the ISMRM abstract
                    % T2 relaxation times of 18 brain metabolites determined in 83 healthy volunteers in vivo
                    % Milan Scheidegger et al. Proc. Intl. Soc. Mag. Reson. Med. 22 (2014)-- MGSaleh
                    N_H_Metab = 2;
                    MM = 1;
                    
                case 'Lac'
                    EditingEfficiency = 0.94; % determined by FID-A simulations (for TE = 140 ms)
                    T1_Metab = 1.50; % Wijnen et al. 2015 (NMR Biomed)
                    T2_Metab = 0.24; % Madan et al. 2015 (MRM) (NB: this was estimated in brain tumors)
                    N_H_Metab = 3;
                    MM = 1;
            end
            
            MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_Quant(ii) = (MRS_struct.out.(vox{kk}).(target{trg}).Area(ii) ./ MRS_struct.out.(vox{kk}).water.Area(ii)) * ...
                (N_H_Water/N_H_Metab) * MM / EditingEfficiency * ...
                (fracGM * concw_GM * (1-exp(-TR/T1w_GM)) * (exp(-TE/T2w_GM)) / ((1-exp(-TR/T1_Metab)) * (exp(-TE/T2_Metab))) + ...
                fracWM * concw_WM * (1-exp(-TR/T1w_WM)) * (exp(-TE/T2w_WM)) / ((1-exp(-TR/T1_Metab)) * (exp(-TE/T2_Metab))) + ...
                fracCSF * concw_CSF * (1-exp(-TR/T1w_CSF)) * (exp(-TE/T2w_CSF)) / ((1-exp(-TR/T1_Metab)) * (exp(-TE/T2_Metab))));
            MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_QuantCorr(ii)         = MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_Quant(ii) / (fracGM + alpha*fracWM);
            MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_QuantNormTissCorr(ii) = MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_Quant(ii) * CorrFactor;
            
        end
    end
    
    % MRS_struct.Quantify.tissalpha = MRS_struct.out.GABAconciu * CorrFactor;
    % probably make a new output - combination of GannetSeg and what did
    % here...
    
    % save .mat
    
    pdfdirname = './GannetQuantify_output';
    if ~exist(pdfdirname,'dir')
        mkdir(pdfdirname)
    end
    
    if MRS_struct.p.mat
        matname = [pdfdirname '/' 'MRS_struct.mat'];
        save(matname,'MRS_struct');
    end
    
    % build output pdf summary
    
    if ishandle(105)
        clf(105) % MM (170831)
    end
    h = figure(105);
    % MM (170831): Open figure in center of screen
    scr_sz = get(0, 'ScreenSize');
    fig_w = 1000;
    fig_h = 707;
    set(h, 'Position', [(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
    set(h,'Color',[1 1 1]);
    figTitle = 'GannetQuantify Output';
    set(gcf,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
    
    % GABA plot
    % a subplot of the voxel on the brain
    % replot of GABA fit spec
    
    % OUTPUTS:
    % data files - .data and .mask
    % GABA CSF corr
    % GABA reln cnsts
    % GABA alpha corr
    % GABA normalize voxel alpha corr
    
    subplot(2,2,4);
    imagesc(squeeze(MRS_struct.mask.img(ii,:,1:round(size(MRS_struct.mask.img,3)/3))));
    colormap('gray');
    caxis([0 0.5])
    axis equal;
    axis tight;
    axis off;
    
    subplot(2,2,1);
    z=abs(MRS_struct.spec.freq-4.1);
    lowerbound=find(min(z)==z);
    z=abs(MRS_struct.spec.freq-2.79);
    upperbound=find(min(z)==z);
    freqbounds=lowerbound:upperbound;
    freq=MRS_struct.spec.freq(freqbounds);
    plot(MRS_struct.spec.freq, real(MRS_struct.spec.(vox{kk}).GABAGlx.diff(ii,:)), ...
                'k', freq, GABAGlxModel(MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:),freq), 'r');
    
    zz=abs(MRS_struct.spec.freq-3.6);
    Glx_right=find(min(zz)==zz);
    %zz=abs(MRS_struct.spec.freq-3.3);
    %GABA_left=find(min(zz)==zz);
    zz=abs(MRS_struct.spec.freq-2.8);
    GABA_right=find(min(zz)==zz);
    %specbaseline = mean(real(SpectraToPlot(1,GABA_right:GABA_left)),2);
    gabaheight = max(abs(real(MRS_struct.spec.(vox{kk}).GABAGlx.diff(ii,Glx_right:GABA_right))),[],2);
    gabaheight = mean(gabaheight);
    
    yaxismax = 2 *gabaheight; % top spec + 2* height of gaba
    yaxismin = -2* gabaheight; % extend 2* gaba heights below zero
    if yaxismax < yaxismin
        dummy=yaxismin;
        yaxismin=yaxismax;
        yaxismax=dummy;
    end
    axis([0 5  yaxismin yaxismax]);
    set(gca,'YTick',[]);
    set(gca,'XLim',[0 4.5]);
    set(gca,'XDir','reverse');
    
    subplot(2,2,2);  % output results
    axis off;
    
    % tmp = ['GABAconc(iu) tissue corr (CSF corrected):  ' num2str(MRS_struct.out.GABAconciuTissCorr(ii))];
    %     text(0, 0.87, tmp, 'HorizontalAlignment', 'left', ...
    %             'VerticalAlignment', 'top',...
    %             'FontName', 'Helvetica','FontSize',13);
    
    tmp = ['GABAconc(iu) with tissue relaxation and vis:  ' num2str(MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_Quant(ii))];
    text(0, 0.77, tmp, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top',...
        'FontName', 'Helvetica','FontSize',13);
    
    tmp = ['GABAconc(iu) alpha-corrected:  ' num2str(MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_QuantCorr(ii))];
    text(0, 0.67, tmp, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top',...
        'FontName', 'Helvetica','FontSize',13);
    
    tmp = '    (uses tissue specific relaxation and visibility constants)';
    text(0, 0.6, tmp, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top',...
        'FontName', 'Helvetica','FontSize',13);
    
    tmp = ['GABAconc(iu) alpha-corrected, average voxel-normalized:  ' num2str(MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_QuantNormTissCorr(ii))];
    text(0, 0.5, tmp, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top',...
        'FontName', 'Helvetica','FontSize',13);
    
    tmp = '    (uses tissue specific relaxation and visibility constants)' ;
    text(0, 0.43, tmp, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top',...
        'FontName', 'Helvetica','FontSize',13);
    
    C = MRS_struct.metabfile{ii};
    if size(C,2) > 30
        [~,y] = fileparts(C);
    else
        y = C;
    end
    tmp = ['MRS data:  ' y];
    tmp = regexprep(tmp, '_','-');
    text(0,0.25, tmp, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top',...
        'FontName', 'Helvetica','FontSize',13);
    
    D = MRS_struct.mask.T1image{ii} ;
    if size(D,2) >30
        [~,y] = fileparts(D);
    else
        y = D;
    end
    tmp = ['Anatomical image:  ' y];
    tmp = regexprep(tmp, '_','-');
    text(0,0.15, tmp, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top',...
        'FontName', 'Helvetica','FontSize',13);
    
    tmp = ['QuantifyVer:  ' MRS_struct.version.quantify];
    text(0,0.0, tmp, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top',...
        'FontName', 'Helvetica','FontSize',13);
    
    subplot(2,2,3);
    rejectframesplot = (1./MRS_struct.out.reject(:,ii).') .*  MRS_struct.fids.waterfreq(ii,:);
    plot(1:size(MRS_struct.fids.data,2), ...
        MRS_struct.fids.waterfreq(ii,:)', '-', ...
        1:size(MRS_struct.fids.data,2), rejectframesplot, 'ro');
    set(gca,'XLim',[0 size(MRS_struct.fids.data,2)]);
    xlabel('time'); ylabel('\omega_0');
    title('Water Frequency, ppm');
    
    subplot(2,2,4);
    axis off;
    script_path=which('GannetFit');
    Gannet_logo=[script_path(1:(end-12)) '/Gannet3_logo.png'];
    A2=imread(Gannet_logo,'png','BackgroundColor',[1 1 1]);
    axes('Position',[0.80, 0.05, 0.15, 0.15]);
    image(A2); axis off; axis square;
    
    %%%% Save PDF %%%%%
    pfil_nopath = MRS_struct.metabfile{ii};
    fullpath = MRS_struct.metabfile{ii};
    
    if(strcmp(fullpath(1:2) , '..'))
        fullpath = fullpath(4:end);
    end
    
    if strcmpi(MRS_struct.p.vendor,'Philips_data')
        fullpath = MRS_struct.metabfile{ii};
        if strcmp(fullpath(1:2) , '..')
            fullpath = fullpath(4:end);
        end
        fullpath = regexprep(fullpath, '.data', '_data');
        fullpath = regexprep(fullpath, '\', '_');
        fullpath = regexprep(fullpath, '/', '_');
    end
    
    tmp = strfind(pfil_nopath,'/');
    tmp2 = strfind(pfil_nopath,'\');
    if(tmp)
        lastslash=tmp(end);
    elseif (tmp2)
        %maybe it's Windows...
        lastslash=tmp2(end);
    else
        % it's in the current dir...
        lastslash=0;
    end
    if(strcmpi(MRS_struct.p.vendor,'Philips'))
        tmp = strfind(pfil_nopath, '.sdat');
        tmp1= strfind(pfil_nopath, '.SDAT');
        if size(tmp,1)>size(tmp1,1)
            dot7 = tmp(end); % just in case there's another .sdat somewhere else...
        else
            dot7 = tmp1(end); % just in case there's another .sdat somewhere else...
        end
    elseif(strcmpi(MRS_struct.p.vendor,'GE'))
        tmp = strfind(pfil_nopath, '.7');
        dot7 = tmp(end); % just in case there's another .7 somewhere else...
    elseif(strcmpi(MRS_struct.p.vendor,'Philips_data'))  % make this be sdat
        tmp = strfind(pfil_nopath, '.data');
        dot7 = tmp(end); % just in case there's another .data somewhere else...
    end
    pfil_nopath = pfil_nopath(lastslash+1:dot7-1);
    
    % Save PDF output
    set(gcf, 'PaperUnits', 'inches');
    set(gcf,'PaperSize',[11 8.5]);
    set(gcf,'PaperPosition',[0 0 11 8.5]);
    if(strcmpi(MRS_struct.p.vendor,'Philips_data'))
        pdfname = [pdfdirname '/' fullpath '_quantify.pdf'];
    else
        pdfname = [pdfdirname '/' pfil_nopath  '_quantify.pdf'];
    end
    saveas(gcf, pdfname);
    
end
