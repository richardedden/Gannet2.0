function MRS_struct = GannetSegment(MRS_struct)

% Relies on SPM being installed
%
% Runs segmentation script if segmented images not present according to
% file convention of c1, c2 and c3 as prefixes on the anatomical image name
% for the GM, WM and CSF segmentations. If these files are present, they
% are loaded and used for the voxel segmentation

MRS_struct.version.segment = '170914';

if MRS_struct.p.PRIAM % deciding how many regions are there -- MGSaleh 2016
    vox = {MRS_struct.p.Vox};
else
    vox = {MRS_struct.p.Vox{1}};
end

% Set up SPM for batch processing
spm('defaults','fmri');
spm_jobman('initcfg');

for ii = 1:length(MRS_struct.metabfile)
    
    % 1 - Take nifti from GannetCoRegister and segment it in SPM
    
    [T1dir, T1name, T1ext] = fileparts(MRS_struct.mask.T1image{ii});
    anatimage = MRS_struct.mask.T1image{ii};
    
    % Check to see if segmentation already done - if not, do it
    % Check which SPM version is installed and segment accordingly
    tmp = [T1dir '/c1' T1name T1ext];
    if ~exist(tmp,'file')
        spmversion = fileparts(which('spm'));
        if strcmpi(spmversion(end-4:end),'spm12')
            CallSPM12segmentation(anatimage);
        else
            CallSPM8segmentation(anatimage);
        end
    end
    
    % 2 - Determine GM, WM and CSF fractions for each voxel
    
    if strcmp(T1dir,'')
        T1dir='.';
    end
    
    GM  = [T1dir '/c1' T1name T1ext];
    WM  = [T1dir '/c2' T1name T1ext];
    CSF = [T1dir '/c3' T1name T1ext];
    
    GMvol  = spm_vol(GM);
    WMvol  = spm_vol(WM);
    CSFvol = spm_vol(CSF);
    
    % Loop over voxels if PRIAM
    for kk = 1:length(vox)
        
        voxmaskvol = spm_vol(cell2mat(MRS_struct.mask.outfile(ii)));
        
        % GM
        O_GMvox.fname = [T1dir '/c1' T1name '_GM_mask.nii'];
        O_GMvox.descrip = 'GMmasked_MRS_Voxel_Mask';
        O_GMvox.dim = voxmaskvol.dim;
        O_GMvox.dt = voxmaskvol.dt;
        O_GMvox.mat = voxmaskvol.mat;
        GM_voxmask_vol = GMvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
        O_GMvox = spm_write_vol(O_GMvox, GM_voxmask_vol);
        
        % WM
        O_WMvox.fname = [T1dir '/c2' T1name '_WM_mask.nii'];
        O_WMvox.descrip = 'WMmasked_MRS_Voxel_Mask';
        O_WMvox.dim = voxmaskvol.dim;
        O_WMvox.dt = voxmaskvol.dt;
        O_WMvox.mat = voxmaskvol.mat;
        WM_voxmask_vol = WMvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
        O_WMvox = spm_write_vol(O_WMvox, WM_voxmask_vol);
        
        % CSF
        O_CSFvox.fname = [T1dir '/c3' T1name '_CSF_mask.nii'];
        O_CSFvox.descrip = 'CSFmasked_MRS_Voxel_Mask';
        O_CSFvox.dim = voxmaskvol.dim;
        O_CSFvox.dt = voxmaskvol.dt;
        O_CSFvox.mat = voxmaskvol.mat;
        CSF_voxmask_vol = CSFvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
        O_CSFvox = spm_write_vol(O_CSFvox, CSF_voxmask_vol);
        
        % ***NB***
        % For a subject with multiple voxels, the segmented voxels will
        % get overwritten
        
        % 3 - Calculate an adjusted gabaiu and output it to the structure
        
        gmsum = sum(sum(sum(O_GMvox.private.dat(:,:,:))));
        wmsum = sum(sum(sum(O_WMvox.private.dat(:,:,:))));
        csfsum = sum(sum(sum(O_CSFvox.private.dat(:,:,:))));
        
        gmfra = gmsum/(gmsum+wmsum+csfsum);
        wmfra = wmsum/(gmsum+wmsum+csfsum);
        csffra = csfsum/(gmsum+wmsum+csfsum);
        
        tissuefra = gmfra+wmfra;
        
        MRS_struct.out.(vox{kk}).tissue.GMfra(ii) = gmfra;
        MRS_struct.out.(vox{kk}).tissue.WMfra(ii) = wmfra;
        MRS_struct.out.(vox{kk}).tissue.CSFfra(ii) = csffra;
        
        % Correction of institutional units only feasible if water-scaling is
        % performed, skip otherwise (GO 07/13/2017)
        % MM (170831): Loop over voxels and metabolites
        if strcmp(MRS_struct.p.Reference_compound,'H2O')
            if MRS_struct.p.HERMES
                target = {MRS_struct.p.target, MRS_struct.p.target2};
            else
                target = {MRS_struct.p.target};
            end
            
            for trg = 1:length(target)
                if strcmp(target{trg},'GABAGlx')
                    MRS_struct.out.(vox{kk}).GABA.ConcIU_CSFcorr(ii) = ...
                        MRS_struct.out.(vox{kk}).GABA.ConcIU(ii) / tissuefra;
                    MRS_struct.out.(vox{kk}).Glx.ConcIU_CSFcorr(ii) = ...
                        MRS_struct.out.(vox{kk}).Glx.ConcIU(ii) / tissuefra;
                else
                    MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_CSFcorr(ii) = ...
                        MRS_struct.out.(vox{kk}).(target{trg}).ConcIU(ii) / tissuefra;
                end
            end
        end
        
        % 4 - Build output
        
        if ishandle(104)
            clf(104); % MM (170831)
        end
        h = figure(104);
        % MM (170831): Open figure in center of screen
        scr_sz = get(0, 'ScreenSize');
        fig_w = 1000;
        fig_h = 707;
        set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
        set(h,'Color',[1 1 1]);
        figTitle = 'GannetSegment Output';
        set(gcf,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
        
        % Voxel co-registration
        subplot(2,2,1);
        size_max = size(MRS_struct.mask.img{ii},1);
        imagesc(MRS_struct.mask.img{ii}(:,size_max+(1:size_max)));
        colormap('gray');
        caxis([0 0.5]);
        axis equal;
        axis tight;
        axis off;
        
        % MM (180112)
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
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
        
        text_pos = 1;
        
        % Print correction of institutional units only feasible if water-scaling is
        % performed, skip otherwise (GO 07/13/2017)
        if strcmp(MRS_struct.p.Reference_compound,'H2O')
            for trg = 1:length(target)
                
                text_pos = text_pos - 0.1;
                
                switch target{trg}
                    case 'GABA'
                        tmp1 = 'GABA+ (CSF-corrected)';
                        tmp2 = sprintf(': %.3g i.u.', MRS_struct.out.(vox{kk}).GABA.ConcIU_CSFcorr(ii));
                        
                    case {'Glx','GSH','Lac'}
                        tmp1 = [target{trg} ' (CSF-corrected)'];
                        tmp2 = sprintf(': %.3g i.u.', MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_CSFcorr(ii));
                        
                    case 'GABAGlx'
                        tmp1 = 'GABA+/Glx (CSF-corrected)';
                        tmp2 = sprintf(': %.3g/%.3g i.u.', MRS_struct.out.(vox{kk}).GABA.ConcIU_CSFcorr(ii), ...
                            MRS_struct.out.(vox{kk}).Glx.ConcIU_CSFcorr(ii));
                end
                
                text(0, text_pos, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
                text(0.5, text_pos, tmp2, 'FontName', 'Helvetica', 'FontSize', 10);
                
            end
        end
        
        tmp1 = 'GM voxel fraction';
        tmp2 = sprintf(': %.3g', MRS_struct.out.(vox{kk}).tissue.GMfra(ii));
        text(0, text_pos-0.1, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
        text(0.5, text_pos-0.1, tmp2, 'FontName', 'Helvetica', 'FontSize', 10);
        
        tmp1 = 'WM voxel fraction';
        tmp2 = sprintf(': %.3g', MRS_struct.out.(vox{kk}).tissue.WMfra(ii));
        text(0, text_pos-0.2, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
        text(0.5, text_pos-0.2, tmp2, 'FontName', 'Helvetica', 'FontSize', 10);
        
        tmp1 = 'CSF voxel fraction';
        tmp2 = sprintf(': %.3g', MRS_struct.out.(vox{kk}).tissue.CSFfra(ii));
        text(0, text_pos-0.3, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
        text(0.5, text_pos-0.3, tmp2, 'FontName', 'Helvetica', 'FontSize', 10);
        
        % MM (180112)
        tmp1 = 'Filename';
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{ii*2-1});
        else
            [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{ii});
        end
        text(0, text_pos-0.4, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
        text(0.5, text_pos-0.4, [': ' tmp2 tmp3], 'FontName', 'Helvetica', 'FontSize', 10, 'Interpreter', 'none');
        
        tmp1 = 'Anatomical image';
        [~,tmp2,tmp3] = fileparts(MRS_struct.mask.T1image{ii}); % MM (180112)
        text(0, text_pos-0.5, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
        text(0.5, text_pos-0.5, [': ' tmp2 tmp3], 'FontName', 'Helvetica', 'FontSize', 10, 'Interpreter', 'none');
        
        tmp1 = 'SegmentVer';
        tmp2 = [': ' MRS_struct.version.segment];
        text(0, text_pos-0.6, tmp1, 'FontName', 'Helvetica', 'FontSize', 10);
        text(0.5, text_pos-0.6, tmp2, 'FontName', 'Helvetica', 'FontSize', 10);
        
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
        if ~exist('GannetSegment_output','dir')
            mkdir GannetSegment_output;
        end
        
        % For Philips .data
        if strcmpi(MRS_struct.p.vendor,'Philips_data')
            fullpath = MRS_struct.metabfile{ii};
            fullpath = regexprep(fullpath, '.data', '_data');
            fullpath = regexprep(fullpath, '\', '_');
            fullpath = regexprep(fullpath, '/', '_');
        end
        
        % MM (180112)
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii*2-1});
        else
            [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii});
        end
        
        % Save PDF output
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperSize',[11 8.5]);
        set(gcf,'PaperPosition',[0 0 11 8.5]);
        if strcmpi(MRS_struct.p.vendor,'Philips_data')
            pdfname = fullfile('GannetSegment_output', [fullpath '_segment.pdf']); % MM (180112)
        else
            pdfname = fullfile('GannetSegment_output', [metabfile_nopath '_segment.pdf']); % MM (180112)
        end        
        saveas(gcf, pdfname);
        
    end
    
end

end



