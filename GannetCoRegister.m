function MRS_struct = GannetCoRegister(MRS_struct, nii_name, rot_folder)

%Coregistration of MRS voxel volumes to imaging datasets, based on headers.

MRS_struct.versioncoreg = '170831';

if MRS_struct.ii ~= length(nii_name)
    error('The number of nifti files does not match the number of MRS files processed by GannetLoad.');
end

for ii = 1:length(MRS_struct.gabafile)
    
    %Ultimately this switch will not be necessary...
    switch MRS_struct.p.vendor
        
        case 'Philips'
            fname = MRS_struct.gabafile{ii};
            sparname = [fname(1:(end-4)) MRS_struct.p.spar_string];
            MRS_struct = GannetMask_Philips(sparname, nii_name{ii}, MRS_struct, ii);
            
        case 'Philips_data'
            if exist(MRS_struct.gabafile_sdat,'file') % MM (170720)
                MRS_struct.p.vendor = 'Philips';
                MRS_struct.gabafile_data = MRS_struct.gabafile;
                MRS_struct.gabafile_data = MRS_struct.gabafile;
                MRS_struct.gabafile = MRS_struct.gabafile_sdat;
                MRS_struct = GannetCoRegister(MRS_struct,nii_name);
                MRS_struct.gabafile = MRS_struct.gabafile_data;
                MRS_struct.p.vendor = 'Philips_data';
            else
                error([MRS_struct.p.vendor ' format does not include voxel location information in the header. See notes in GannetCoRegister.']);
                %If this comes up, once GannetLoad has been read:
                %1. Switch vendor to Philips
                %       MRS_struct.p.vendor = 'Philips';
                %2. Copy .data filenames.
                %       MRS_struct.gabafile_data = MRS_struct.gabafile;
                %3. Replace the list with the corrsponding SDAT files (in correct order)
                %        MRS_struct.gabafile = {'SDATfile1.sdat' 'SDATfile2.SDAT'};
                %4. Rerun GannetCoRegister
                %
                %5.  Copy .sdat filenames and replace .data ones. Tidy up.
                %       MRS_struct.gabafile_sdat = MRS_struct.gabafile;
                %       MRS_struct.gabafile = MRS_struct.gabafile_data;
                %       MRS_struct.p.vendor = 'Philips_data'
            end
            
        case 'Siemens'
            fname = MRS_struct.gabafile{ii};
            MRS_struct = GannetMask_Siemens(fname, nii_name{ii}, MRS_struct, ii);
            
        case 'Siemens_twix'
            fname = MRS_struct.gabafile{ii};
            MRS_struct = GannetMask_SiemensTWIX(fname, nii_name{ii}, MRS_struct, ii);
            
        case 'GE'
            fname = MRS_struct.gabafile{ii};
            if ~exist('rot_folder','var')
                rot_folder = nii_name;
            end
            MRS_struct = GannetMask_GE(fname, nii_name{ii}, MRS_struct, rot_folder{ii}, ii);
            
    end
    
    % Build output figure
    if ishandle(103)
        clf(103) % MM (170720)
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
    text(0.5,0.75, tmp,'HorizontalAlignment','right', ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    D = MRS_struct.mask.outfile{ii};
    if size(D,2) > 30
        [~,tmp,~] = fileparts(D);
    else
        tmp = D;
    end
    tmp = [': ' regexprep(tmp,'_','-')];
    text(0.5,0.75, tmp, ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    
    tmp = 'Spatial parameters ';
    text(0.5,0.63, tmp, 'HorizontalAlignment','right', ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    tmp = ': [LR, AP, FH]';
    text(0.5, 0.63, tmp, ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    
    tmp = 'Dimension ';
    text(.5,0.51, tmp, 'HorizontalAlignment','right', ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    tmp = [': [' num2str(MRS_struct.p.voxdim(ii,1)) ', ' num2str(MRS_struct.p.voxdim(ii,2)) ', ' num2str(MRS_struct.p.voxdim(ii,3)) '] mm'];
    text(.5,0.51, tmp, ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    
    tmp = 'Volume ';
    text(.5,0.39, tmp, 'HorizontalAlignment','right', ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    vol = MRS_struct.p.voxdim(ii,1)*MRS_struct.p.voxdim(ii,2)*MRS_struct.p.voxdim(ii,3)*.001;
    tmp = [': ' num2str(vol) ' mL'];
    text(.5, 0.39, tmp, ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    
    tmp = 'Position ';
    text(0.5,0.27, tmp, 'HorizontalAlignment','right', ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    tmp = [': [' num2str(MRS_struct.p.voxoff(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxoff(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxoff(ii,3), '%3.1f') '] mm'];
    text(.5, 0.27, tmp, ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    
    tmp = 'Angulation ';
    text(0.5,0.15, tmp, 'HorizontalAlignment','right', ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    tmp = [': [' num2str(MRS_struct.p.voxang(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,3), '%3.1f') '] deg'];
    text(.5, 0.15, tmp, ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    
    tmp = 'CoRegVer ';
    text(0.5,0.03, tmp, 'HorizontalAlignment','right', ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    tmp = [': ' MRS_struct.versioncoreg];
    text(.5, 0.03, tmp, ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica','FontSize',13);
    
    h=subplot(2,3,1:3);
    
    t = ['Voxel from ' MRS_struct.gabafile{ii} ' on ' MRS_struct.mask.T1image(ii,:)];
    t = regexprep(t, '_','-');
    
    imagesc(squeeze(MRS_struct.mask.img(ii,:,:)));
    colormap('gray');
    caxis([0 .5]) % range of 0 to 0.5 seems to work best for now - could calc optimal range later
    axis equal;
    axis tight;
    axis off;
    text(10,size(MRS_struct.mask.img,2)/2,'L','Color',[1 1 1]);
    text(size(MRS_struct.mask.img,3)-15,size(MRS_struct.mask.img,2)/2,'R','Color',[1 1 1]);
    get(h,'pos'); % get position of axes
    set(h,'pos',[0.0 0.15 1 1]) % move the axes slightly
    title(t, 'FontName', 'Helvetica','FontSize',15);
    
    script_path=which('GannetLoad');
    Gannet_logo=[script_path(1:(end-13)) '/Gannet3_logo.png'];
    A2=imread(Gannet_logo,'png','BackgroundColor',[1 1 1]);
    axes('Position',[0.80, 0.08, 0.15, 0.15]);
    image(A2);axis off; axis square;
    
    %%%% Save PDF %%%%%
    pdfdirname = 'GannetCoRegister_output'; % MM (170720)
    pfil_nopath = MRS_struct.gabafile{ii};
    
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
    elseif(strcmpi(MRS_struct.p.vendor,'Philips_data'))
        tmp = strfind(pfil_nopath, '.data');
        dot7 = tmp(end); % just in case there's another .data somewhere else...
    elseif(strcmpi(MRS_struct.p.vendor,'Siemens'))
        tmp = strfind(pfil_nopath, '.rda');
        dot7 = tmp(end); % just in case there's another .rda somewhere else...
    elseif(strcmpi(MRS_struct.p.vendor,'Siemens_twix'))
        tmp = strfind(pfil_nopath, '.dat');
        dot7 = tmp(end); % just in case there's another .dat somewhere else...
    end
    
    pfil_nopath = pfil_nopath(lastslash+1:dot7-1);
    
    % Save PDF output
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize',[11 8.5]);
    set(gcf,'PaperPosition',[0 0 11 8.5]);
    if strcmpi(MRS_struct.p.vendor,'Philips_data')
        pdfname = [pdfdirname '/' fullpath '_coreg.pdf'];
    else
        pdfname = [pdfdirname '/' pfil_nopath  '_coreg.pdf'];
    end
    if ~exist(pdfdirname,'dir')
        mkdir(pdfdirname)
    end
    saveas(gcf, pdfname);    
    
end


