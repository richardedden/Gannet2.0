function MRS_struct = GannetQuantify(MRS_struct)

MRS_struct.version.quantify = '170915';

cWM = 1; % concentration of GABA in pure WM
cGM = 2; % concentration of GABA in pure GM
% NB: cWM/cGM is term used ratio not absolute value is important piece
alpha = cWM/cGM;

% Constants
% From Wansapura et al. 1999 (JMRI)
%        T1          T2
% WM   832 +/- 10  79.2 +/- 0.6
% GM  1331 +/- 13  110 +/- 2
%
% From Lu et al. 2005 (JMRI)
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
% Gasparovic et al. 2006 (MRM) uses relative densities, ref to
% Ernst et al. 1993 (JMR)
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

for ii = 1:length(MRS_struct.metabfile)
    
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
    
    TR = MRS_struct.p.TR(ii)/1000;
    TE = MRS_struct.p.TE(ii)/1000;
    
    for kk = 1:length(vox)
        
        meanGMfra = mean(MRS_struct.out.(vox{kk}).tissue.GMfra); % average GM fraction across subjects
        meanWMfra = mean(MRS_struct.out.(vox{kk}).tissue.WMfra); % average WM fraction across subjects
        
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
            
            MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_TissCorr(ii) = ...
                (MRS_struct.out.(vox{kk}).(target{trg}).Area(ii) ./ MRS_struct.out.(vox{kk}).water.Area(ii)) * ...
                (N_H_Water/N_H_Metab) * MM / EditingEfficiency * ...
                (fracGM * concw_GM * (1-exp(-TR/T1w_GM)) * (exp(-TE/T2w_GM)) / ((1-exp(-TR/T1_Metab)) * (exp(-TE/T2_Metab))) + ...
                fracWM * concw_WM * (1-exp(-TR/T1w_WM)) * (exp(-TE/T2w_WM)) / ((1-exp(-TR/T1_Metab)) * (exp(-TE/T2_Metab))) + ...
                fracCSF * concw_CSF * (1-exp(-TR/T1w_CSF)) * (exp(-TE/T2w_CSF)) / ((1-exp(-TR/T1_Metab)) * (exp(-TE/T2_Metab))));
            MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_AlphaTissCorr(ii) = MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_TissCorr(ii) / (fracGM + alpha*fracWM);
            MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_AlphaTissCorr_GrpNorm(ii) = MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_TissCorr(ii) * CorrFactor;
            
        end
    end
    
    % Build output figure    
    if ishandle(105)
        clf(105); % MM (170831)
    end
    h = figure(105);
    % MM (170831): Open figure in center of screen
    scr_sz = get(0, 'ScreenSize');
    fig_w = 1000;
    fig_h = 707;
    set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
    set(h,'Color',[1 1 1]);
    figTitle = 'GannetQuantify Output';
    set(gcf,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
    
    % Voxel co-registration
    subplot(2,2,1);
    imagesc(squeeze(MRS_struct.mask.img(ii,:,1:round(size(MRS_struct.mask.img,3)/3))));
    colormap('gray');
    caxis([0 0.5]);
    axis equal;
    axis tight;
    axis off;
    
    % MM (180112)
    if strcmp(MRS_struct.p.vendor,'Siemens')
        [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii*2-1});
    else
        [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii});
    end
    [~,tmp3,tmp4] = fileparts(MRS_struct.mask.T1image{ii});
    t = ['Voxel from ' tmp tmp2 ' on ' tmp3 tmp4];
    title(t, 'Interpreter', 'none');
    
    % Post-alignment spectra + model fits
    subplot(2,2,3);
    GannetPlotPrePostAlign2(MRS_struct, vox, ii);
    if MRS_struct.p.HERMES
        title('Edited Spectra and Model Fits');
    else
        title('Edited Spectrum and Model Fit');
    end
    xlabel('ppm');
    set(gca,'YTick',[]);
    
    % Output results
    subplot(2,2,2);
    axis off;
    
    if MRS_struct.p.HERMES
        target = {MRS_struct.p.target, MRS_struct.p.target2};
    else
        target = {MRS_struct.p.target};
    end
    
    for trg = 1:length(target)
        
        switch target{trg}
            case 'GABA'
                tmp2 = 'GABA+';
            case {'Glx','GSH','Lac'}
                tmp2 = target{trg};
            case 'GABAGlx'
                tmp2 = 'GABA+/Glx';
        end
        
        shift = 0;
        
        for jj = 1:3
            
            text_pos = 0.9;
            
            if jj == 1
                tmp1 = 'Relaxation-, tissue-corrected:';
                if strcmp(target{trg},'GABAGlx')
                    tmp3 = sprintf(': %.3g/%.3g i.u.', MRS_struct.out.(vox{kk}).GABA.ConcIU_TissCorr(ii), ...
                        MRS_struct.out.(vox{kk}).Glx.ConcIU_TissCorr(ii));
                else
                    tmp3 = sprintf(': %.3g i.u.', MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_TissCorr(ii));
                end
            elseif jj == 2
                text_pos = text_pos - 0.2 - shift;
                tmp1 = 'Relaxation-, tissue-, alpha-corrected:';
                if strcmp(target{trg},'GABAGlx')
                    tmp3 = sprintf(': %.3g/%.3g i.u.', MRS_struct.out.(vox{kk}).GABA.ConcIU_AlphaTissCorr(ii), ...
                        MRS_struct.out.(vox{kk}).Glx.ConcIU_AlphaTissCorr(ii));
                else
                    tmp3 = sprintf(': %.3g i.u.', MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_AlphaTissCorr(ii));
                end
            elseif jj == 3
                text_pos = text_pos - 0.4 - shift;
                tmp1 = 'Relaxation-, tissue-, alpha-corrected (average-voxel-normalized):';
                if strcmp(target{trg},'GABAGlx')
                    tmp3 = sprintf(': %.3g/%.3g i.u.', MRS_struct.out.(vox{kk}).GABA.ConcIU_AlphaTissCorr_GrpNorm(ii), ...
                        MRS_struct.out.(vox{kk}).Glx.ConcIU_AlphaTissCorr_GrpNorm(ii));
                else
                    tmp3 = sprintf(': %.3g i.u.', MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_AlphaTissCorr_GrpNorm(ii));
                end
            end
            
            text(0, text_pos, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
            text_pos = text_pos - 0.1*trg;
            text(0, text_pos, tmp2, 'FontName', 'Helvetica', 'FontSize', 10);
            text(0.375, text_pos, tmp3, 'FontName', 'Helvetica', 'FontSize', 10);
            
            if MRS_struct.p.HERMES
                shift = shift + 0.1*(numel(target)-1);
            end
            
        end
    end
    
    % MM (180112)
    tmp1 = 'Filename';
    if strcmp(MRS_struct.p.vendor,'Siemens')
        [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{ii*2-1});
    else
        [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{ii});
    end
    text(0, text_pos-0.1, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
    text(0.375, text_pos-0.1, [': ' tmp2 tmp3], 'FontName', 'Helvetica', 'FontSize', 10, 'Interpreter', 'none');
    
    tmp1 = 'Anatomical image';
    [~,tmp2,tmp3] = fileparts(MRS_struct.mask.T1image{ii}); % MM (180112)
    text(0, text_pos-0.2, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
    text(0.375, text_pos-0.2, [': ' tmp2 tmp3], 'FontName', 'Helvetica', 'FontSize', 10, 'Interpreter', 'none');
    
    tmp1 = 'QuantifyVer';
    tmp2 = [': ' MRS_struct.version.quantify];
    text(0, text_pos-0.3, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
    text(0.375, text_pos-0.3, tmp2, 'FontName', 'Helvetica', 'FontSize', 10);
    
    % Gannet logo
    subplot(2,2,4);
    axis off;
    script_path = which('GannetFit');
    Gannet_logo = [script_path(1:(end-12)) '/Gannet3_logo.png'];
    A2 = imread(Gannet_logo,'png','BackgroundColor',[1 1 1]);
    axes('Position',[0.80, 0.05, 0.15, 0.15]);
    image(A2);
    axis off;
    axis square;
    
    % Create output folder
    if ~exist('GannetQuantify_output','dir')
        mkdir GannetQuantify_output;
    end
    
    if MRS_struct.p.mat
        matname = fullfile('GannetQuantify_output','MRS_struct.mat');
        save(matname,'MRS_struct');
    end
    
    % For Philips .data
    if strcmpi(MRS_struct.p.vendor,'Philips_data')
        fullpath = MRS_struct.metabfile{ii};
        fullpath = regexprep(fullpath, '.data', '_data');
        fullpath = regexprep(fullpath, '\', '_');
        fullpath = regexprep(fullpath, '/', '_');
    end
    
    % MM (180112)
    if strcmp(MRS_struct.p.vendor,'Siemens')
        [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii*2-1});
    else
        [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii});
    end
    
    % Save PDF output
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize',[11 8.5]);
    set(gcf,'PaperPosition',[0 0 11 8.5]);
    if strcmpi(MRS_struct.p.vendor,'Philips_data')
        pdfname = fullfile('GannetQuantify_output', [fullpath '_quantify.pdf']); % MM (180112)
    else
        pdfname = fullfile('GannetQuantify_output', [metabfile_nopath '_quantify.pdf']); % MM (180112)
    end
    saveas(gcf, pdfname);
    
end

