
function MRS_struct = GannetCoRegister(MRS_struct,nii_name,rot_folder)

%Coregistration of MRS voxel volumes to imaging datasets, based on headers. 

MRS_struct.p.coreg = 1;

if (MRS_struct.ii ~= length(nii_name))
       error('The number of nifti files does not match the number of MRS files processed by GannetLoad.'); 
end

for ii = 1:MRS_struct.ii
    ii
    nii_name{ii} 

%Ultimately this switch will not be necessary...
    switch MRS_struct.p.vendor
    
    case 'Philips'
            fname = MRS_struct.gabafile{ii};
            sparname = [fname(1:(end-4)) MRS_struct.p.spar_string];
            MRS_struct=GannetMask_Philips(sparname,nii_name{ii},MRS_struct, ii);
    case 'Philips_data'
        if exist(MRS_struct.gabafile_sdat)
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
        error(['GannetCoRegister does not yet support ' MRS_struct.p.vendor ' data.']);        
    case 'Siemens_twix'
        error(['GannetCoRegister does not yet support ' MRS_struct.p.vendor ' data.']);        
    case 'GE'
            fname = MRS_struct.gabafile{ii};
            %sparname = [fname(1:(end-4)) MRS_struct.p.spar_string];
            MRS_struct=GannetMask_GE(fname,nii_name{ii},MRS_struct,rot_folder);
    end
    
    
    
    %Build output figure
    h=figure(103);
        set(h, 'Position', [100, 100, 1000, 707]);
        set(h,'Color',[1 1 1]);
        figTitle = ['GannetCoRegister Output'];
        set(gcf,'Name',figTitle,'Tag',figTitle, 'NumberTitle','off');
              

        h=subplot(2,2,1:2);
        size(squeeze(MRS_struct.mask.img(MRS_struct.ii,:,:)))
        imagesc(squeeze(MRS_struct.mask.img(MRS_struct.ii,:,:)));
        colormap('gray');
        caxis([0 1])
        axis equal;
        axis tight;
        axis off;
        text(10,size(MRS_struct.mask.img,2)/2,'L','Color',[1 1 1]);
        text(size(MRS_struct.mask.img,3)-15,size(MRS_struct.mask.img,2)/2,'R','Color',[1 1 1]);
        p = get(h,'pos'); % get position of axes
        set(h,'pos',[0.05 0.4 0.9 0.45]) % move the axes slightly
        
        script_path=which('GannetLoad');
              % CJE update for GE
    %          Gannet_circle=[script_path(1:(end-12)) 'GANNET_circle.png'];
              Gannet_circle_white=[script_path(1:(end-13)) '/GANNET_circle_white.jpg'];
    %          A=imread(Gannet_circle);
              A2=imread(Gannet_circle_white);
              hax=axes('Position',[0.80, 0.08, 0.15, 0.15]);
              %set(gca,'Units','normalized');set(gca,'Position',[0.05 0.05 1.85 0.15]);
              image(A2);axis off; axis square;

        %axis off;

        
        
           %%%%  Save EPS %%%%%
           
           %110624
epsdirname = [ './MRSCoReg_' datestr(clock,'yymmdd') ];

           
%     if strcmp(MRS_struct.p.vendor,'Siemens')
%     pfil_nopath = MRS_struct.gabafile{ii*2-1};
%     else

pfil_nopath = MRS_struct.gabafile{ii};
%     end

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
%     elseif(strcmpi(MRS_struct.p.vendor,'Siemens'))
%         tmp = strfind(pfil_nopath, '.rda');
%         dot7 = tmp(end); % just in case there's another .data somewhere else...
%     elseif(strcmpi(MRS_struct.p.vendor,'Siemens_twix'))
%         tmp = strfind(pfil_nopath, '.dat');
%         dot7 = tmp(end); % just in case there's another .dat somewhere else...
    end
    pfil_nopath = pfil_nopath( (lastslash+1) : (dot7-1) );
    %Save pdf output
    set(gcf, 'PaperUnits', 'inches');
    set(gcf,'PaperSize',[11 8.5]);
    set(gcf,'PaperPosition',[0 0 11 8.5]);
    if(strcmpi(MRS_struct.p.vendor,'Philips_data'))
        pdfname=[ epsdirname '/' fullpath '.pdf' ];
    else
        pdfname=[ epsdirname '/' pfil_nopath  '.pdf' ];
    end
    %epsdirname
    if(exist(epsdirname,'dir') ~= 7)
        epsdirname
        mkdir(epsdirname)
    end
    saveas(gcf, pdfname);
        
        
end