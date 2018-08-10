function MRS_struct = GannetSegment(MRS_struct)

% Relies on SPM12 being installed
%
% Runs segmentation script if segmented images not present according to
% file convention of c1, c2 and c3 as prefixes on the anatomical image name
% for the GM, WM and CSF segmentations. If these files are present, they
% are loaded and used for the voxel segmentation

MRS_struct.version.segment = '180807';

% First check if SPM12 is installed and on the search path
spmversion = fileparts(which('spm'));
if isempty(spmversion)
    error('SPM not found! Please install SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12) and make sure it is on your search path.');
elseif strcmpi(spmversion(end-3:end),'spm8')
    error(['SPM8 detected! Gannet 3.0 no longer supports SPM8. ' ...
           'Please install SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12) and make sure it is on your search path.']);
end

if MRS_struct.p.PRIAM % deciding how many regions are there -- MGSaleh 2016
    vox = MRS_struct.p.Vox;
else
    vox = {MRS_struct.p.Vox{1}};
end

numscans = numel(MRS_struct.metabfile);
if strcmpi(MRS_struct.p.vendor,'Siemens_rda')
    numscans = numscans/2;
end

% Set up SPM for batch processing
spm('defaults','fmri');
spm_jobman('initcfg');

for ii = 1:numscans
    
    % Loop over voxels if PRIAM
    for kk = 1:length(vox)
        
        % 1 - Take nifti from GannetCoRegister and segment it in SPM
        
        [T1dir, T1name, T1ext] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        anatimage = MRS_struct.mask.(vox{kk}).T1image{ii};
        
        % Check to see if segmentation already done - if not, do it
        tmp = [T1dir '/c1' T1name T1ext];
        if ~exist(tmp,'file')
            CallSPM12segmentation(anatimage);
        end
        
        % 2 - Determine GM, WM and CSF fractions for each voxel
        
        if strcmp(T1dir,'')
            T1dir = '.';
        end
        
        GM  = [T1dir '/c1' T1name T1ext];
        WM  = [T1dir '/c2' T1name T1ext];
        CSF = [T1dir '/c3' T1name T1ext];
        
        GMvol  = spm_vol(GM);
        WMvol  = spm_vol(WM);
        CSFvol = spm_vol(CSF);
        
        voxmaskvol = spm_vol(cell2mat(MRS_struct.mask.(vox{kk}).outfile(ii)));
        
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
        
        % Output results
        subplot(2,3,4:6);
        axis off;
        
        text_pos = 0.87;
        
        % Print correction of institutional units only feasible if water-scaling is
        % performed, skip otherwise (GO 07/13/2017)
        if strcmp(MRS_struct.p.Reference_compound,'H2O')
            for trg = 1:length(target)
                
                text_pos = text_pos - 0.12;
                
                switch target{trg}
                    case 'GABA'
                        tmp1 = 'GABA+ (CSF-corrected): ';
                        tmp2 = sprintf(' %.3g i.u.', MRS_struct.out.(vox{kk}).GABA.ConcIU_CSFcorr(ii));
                        
                    case {'Glx','GSH','Lac'}
                        tmp1 = [target{trg} ' (CSF-corrected): '];
                        tmp2 = sprintf(': %.3g i.u.', MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_CSFcorr(ii));
                        
                    case 'GABAGlx'
                        tmp1 = 'GABA+/Glx (CSF-corrected): ';
                        tmp2 = sprintf(' %.3g/%.3g i.u.', MRS_struct.out.(vox{kk}).GABA.ConcIU_CSFcorr(ii), ...
                            MRS_struct.out.(vox{kk}).Glx.ConcIU_CSFcorr(ii));
                end
                
                text(0.5, text_pos, tmp1, 'FontName', 'Helvetica', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
                text(0.5, text_pos, tmp2, 'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'FontSize', 13);
                
            end
        end
        
        tmp1 = 'GM voxel fraction: ';
        tmp2 = sprintf(' %.3g', MRS_struct.out.(vox{kk}).tissue.GMfra(ii));
        text(0.5, text_pos-0.12, tmp1, 'FontName', 'Helvetica', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.12, tmp2, 'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        tmp1 = 'WM voxel fraction: ';
        tmp2 = sprintf(' %.3g', MRS_struct.out.(vox{kk}).tissue.WMfra(ii));
        text(0.5, text_pos-0.24, tmp1, 'FontName', 'Helvetica', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.24, tmp2, 'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        tmp1 = 'CSF voxel fraction: ';
        tmp2 = sprintf(' %.3g', MRS_struct.out.(vox{kk}).tissue.CSFfra(ii));
        text(0.5, text_pos-0.36, tmp1, 'FontName', 'Helvetica', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.36, tmp2, 'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        tmp1 = 'Filename: ';
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{ii*2-1});
        else
            [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{ii});
        end
        text(0.5, text_pos-0.48, tmp1, 'FontName', 'Helvetica', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.48, [' ' tmp2 tmp3],  'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'FontSize', 13, 'Interpreter', 'none');
        
        tmp1 = 'Anatomical image: ';
        [~,tmp2,tmp3] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        text(0.5, text_pos-0.6, tmp1, 'FontName', 'Helvetica', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.6, [' ' tmp2 tmp3],  'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'FontSize', 13, 'Interpreter', 'none');
        
        tmp1 = 'SegmentVer: ';
        tmp2 = [' ' MRS_struct.version.segment];
        text(0.5, text_pos-0.72, tmp1, 'FontName', 'Helvetica', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.72, tmp2,  'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        % Voxel segmentation (MM: 180807)
        T1 = spm_read_vols(spm_vol(anatimage));
        img_t     = flipud(voxel2world_space(spm_vol(anatimage), MRS_struct.p.voxoff(ii,:)));
        vox_t     = flipud(voxel2world_space(voxmaskvol, MRS_struct.p.voxoff(ii,:)));
        vox_t_GM  = flipud(voxel2world_space(O_GMvox, MRS_struct.p.voxoff(ii,:)));
        vox_t_WM  = flipud(voxel2world_space(O_WMvox, MRS_struct.p.voxoff(ii,:)));
        vox_t_CSF = flipud(voxel2world_space(O_CSFvox, MRS_struct.p.voxoff(ii,:)));
        img_t = img_t/max(T1(:));
        img_montage = [img_t+0.175*vox_t, img_t+0.21*vox_t_GM, img_t+0.25*vox_t_WM, img_t+0.4*vox_t_CSF];
        
        h = subplot(2,3,1:3);
        imagesc(img_montage);
        colormap('gray');
        img = MRS_struct.mask.(vox{kk}).img{ii};
        img = img(:);
        caxis([0 mean(img(img>0.01)) + 3*std(img(img>0.01))]);
        axis equal;
        axis tight;
        axis off;
        text(floor(size(vox_t,2)/2), 20, 'Voxel', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center'); 
        text(floor(size(vox_t,2)) + floor(size(vox_t,2)/2), 20, 'GM', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center'); 
        text(2*floor(size(vox_t,2)) + floor(size(vox_t,2)/2), 20, 'WM', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center'); 
        text(3*floor(size(vox_t,2)) + floor(size(vox_t,2)/2), 20, 'CSF', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center'); 
        get(h,'pos');
        set(h,'pos',[0 0.15 1 1]);
        
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii*2-1});
        else
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii});
        end
        [~,tmp3,tmp4] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        t = ['Voxel from ' tmp tmp2 ' on ' tmp3 tmp4];
        title(t, 'FontName', 'Helvetica', 'FontSize', 15, 'Interpreter', 'none');
        
        % Gannet logo
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
            pdfname = fullfile('GannetSegment_output', [fullpath '_' vox{kk} '_segment.pdf']); % MM (180112)
        else
            pdfname = fullfile('GannetSegment_output', [metabfile_nopath '_' vox{kk} '_segment.pdf']); % MM (180112)
        end
        saveas(gcf, pdfname);
        
        if ii == numscans
            if MRS_struct.p.mat % save MRS_struct as mat file
                mat_name = ['GannetSegment_output/MRS_struct_' vox{kk} '.mat'];
                save(mat_name,'MRS_struct');
            end
            if MRS_struct.p.csv % export MRS_struct fields into csv file
                if ~strcmp(MRS_struct.p.Reference_compound,'H2O') % this bit of code needed for ExportToCSV when no water ref is provided
                    if MRS_struct.p.HERMES
                        target = {MRS_struct.p.target, MRS_struct.p.target2};
                    else
                        target = {MRS_struct.p.target};
                    end
                end
                ExportToCSV(MRS_struct, target, kk, 'segment');
            end
        end
        
    end
    
end

end



