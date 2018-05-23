function MRS_struct = GannetCoRegister(MRS_struct, nii_name)

%Coregistration of MRS voxel volumes to imaging datasets, based on headers.

MRS_struct.version.coreg = '180523';

if MRS_struct.p.PRIAM % deciding how many voxels there are -- MGSaleh 2016
    vox = MRS_struct.p.Vox;
else
    vox = {MRS_struct.p.Vox{1}};
end

if MRS_struct.ii ~= length(nii_name)
    error('The number of nifti files does not match the number of MRS files processed by GannetLoad.');
end

numscans = numel(MRS_struct.metabfile);
if strcmpi(MRS_struct.p.vendor,'Siemens_rda')
    numscans = numscans/2;
end

for ii = 1:numscans
    % Loop over voxels if PRIAM
    for kk = 1:length(vox)

        %Ultimately this switch will not be necessary...
        switch MRS_struct.p.vendor

            case 'Philips'
                fname = MRS_struct.metabfile{ii};
                sparname = [fname(1:(end-4)) MRS_struct.p.spar_string];
                MRS_struct = GannetMask_Philips(sparname, nii_name{ii}, MRS_struct, ii, vox, kk);

            case 'Philips_data'
                if exist(MRS_struct.metabfile_sdat,'file') % MM (170720)
                    MRS_struct.p.vendor = 'Philips';
                    MRS_struct.metabfile_data = MRS_struct.metabfile;
                    MRS_struct.metabfile = MRS_struct.metabfile_sdat;
                    MRS_struct = GannetCoRegister(MRS_struct,nii_name);
                    MRS_struct.metabfile = MRS_struct.metabfile_data;
                    MRS_struct.p.vendor = 'Philips_data';
                else
                    error([MRS_struct.p.vendor ' format does not include voxel location information in the header. See notes in GannetCoRegister.']);
                    %If this comes up, once GannetLoad has been read:
                    %1. Switch vendor to Philips
                    %       MRS_struct.p.vendor = 'Philips';
                    %2. Copy .data filenames.
                    %       MRS_struct.metabfile_data = MRS_struct.metabfile;
                    %3. Replace the list with the corrsponding SDAT files (in correct order)
                    %        MRS_struct.metabfile = {'SDATfile1.sdat' 'SDATfile2.SDAT'};
                    %4. Rerun GannetCoRegister
                    %
                    %5.  Copy .sdat filenames and replace .data ones. Tidy up.
                    %       MRS_struct.metabfile_sdat = MRS_struct.metabfile;
                    %       MRS_struct.metabfile = MRS_struct.metabfile_data;
                    %       MRS_struct.p.vendor = 'Philips_data'
                end

            case 'Siemens_rda'
                fname = MRS_struct.metabfile{ii*2-1};
                MRS_struct = GannetMask_SiemensRDA(fname, nii_name{ii}, MRS_struct, ii, vox, kk);

            case {'Siemens_twix', 'Siemens_dicom'} 
                fname = MRS_struct.metabfile{ii};
                MRS_struct = GannetMask_SiemensTWIX(fname, nii_name{ii}, MRS_struct, ii, vox, kk);

            case 'GE'
                fname = MRS_struct.metabfile{ii};
                MRS_struct = GannetMask_GE(fname, nii_name{ii}, MRS_struct, ii, vox, kk);

        end

        % Build output figure
        if ishandle(103)
            clf(103); % MM (170720)
        end
        h = figure(103);
        % MM (170629): Open figure in center of screen
        scr_sz = get(0, 'ScreenSize');
        fig_w = 1000;
        fig_h = 707;
        set(h, 'Position', [(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
        set(h,'Color',[1 1 1]);
        figTitle = 'GannetCoRegister Output';
        set(gcf,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');

        subplot(2,3,4:6)
        axis off;

        tmp = 'Mask output ';
        text(0.5, 0.75, tmp,'HorizontalAlignment','right', ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);

        [~,tmp,tmp2] = fileparts(MRS_struct.mask.(vox{kk}).outfile{ii});
        text(0.5, 0.75, [': ' tmp tmp2], ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13,'Interpreter','none');

        tmp = 'Spatial parameters ';
        text(0.5, 0.63, tmp, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);
        tmp = ': [LR, AP, FH]';
        text(0.5, 0.63, tmp, ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);

        tmp = 'Dimension ';
        text(0.5, 0.51, tmp, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);
        tmp = [': [' num2str(MRS_struct.p.voxdim(ii,1)) ', ' num2str(MRS_struct.p.voxdim(ii,2)) ', ' num2str(MRS_struct.p.voxdim(ii,3)) '] mm'];
        text(0.5, 0.51, tmp, ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);

        tmp = 'Volume ';
        text(0.5, 0.39, tmp, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);
        vol = MRS_struct.p.voxdim(ii,1)*MRS_struct.p.voxdim(ii,2)*MRS_struct.p.voxdim(ii,3)*.001;
        tmp = [': ' num2str(vol) ' mL'];
        text(0.5, 0.39, tmp, ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);

        tmp = 'Position ';
        text(0.5, 0.27, tmp, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);
        tmp = [': [' num2str(MRS_struct.p.voxoff(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxoff(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxoff(ii,3), '%3.1f') '] mm'];
        text(0.5, 0.27, tmp, ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);

        tmp = 'Angulation ';
        text(0.5, 0.15, tmp, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);
        tmp = [': [' num2str(MRS_struct.p.voxang(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,3), '%3.1f') '] deg'];
        text(0.5, 0.15, tmp, ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);

        tmp = 'CoRegVer ';
        text(0.5, 0.03, tmp, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);
        tmp = [': ' MRS_struct.version.coreg];
        text(0.5, 0.03, tmp, ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Helvetica','FontSize',13);

        h = subplot(2,3,1:3);

        % MM (180112)
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii*2-1});
        else
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii});
        end
        [~,tmp3,tmp4] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        t = ['Voxel from ' tmp tmp2 ' on ' tmp3 tmp4];

        imagesc(squeeze(MRS_struct.mask.(vox{kk}).img{ii}));
        colormap('gray');
        caxis([0 0.5]); % range of 0 to 0.5 seems to work best for now - could calc optimal range later
        axis equal;
        axis tight;
        axis off;
        text(10,size(MRS_struct.mask.(vox{kk}).img{ii},1)/2,'L','Color',[1 1 1],'FontSize',20);
        text(size(MRS_struct.mask.(vox{kk}).img{ii},2)-20,size(MRS_struct.mask.(vox{kk}).img{ii},1)/2,'R','Color',[1 1 1],'FontSize',20);
        get(h,'pos'); % get position of axes
        set(h,'pos',[0.0 0.15 1 1]) % move the axes slightly
        title(t, 'FontName', 'Helvetica', 'FontSize', 15, 'Interpreter', 'none');

        script_path = which('GannetLoad');
        Gannet_logo = [script_path(1:(end-13)) '/Gannet3_logo.png'];
        A2 = imread(Gannet_logo,'png','BackgroundColor',[1 1 1]);
        axes('Position',[0.80, 0.08, 0.15, 0.15]);
        image(A2);
        axis off;
        axis square;

        % Create output folder
        if ~exist('GannetCoRegister_output','dir')
            mkdir GannetCoRegister_output;
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
            pdfname = fullfile('GannetCoRegister_output', [fullpath '_' vox{kk} '_coreg.pdf']); % MM (180112)
        else
            pdfname = fullfile('GannetCoRegister_output', [metabfile_nopath '_' vox{kk} '_coreg.pdf']); % MM (180112)
        end
        saveas(gcf, pdfname);
        
        % Save MRS_struct as mat file
        if ii == numscans && MRS_struct.p.mat
            % Set up filename
            mat_name = ['GannetCoRegister_output/MRS_struct_' vox{kk} '.mat'];
            save(mat_name,'MRS_struct');
        end
        
    end
    
end




