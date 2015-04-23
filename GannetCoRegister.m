
function MRS_struct = GannetCoRegister(MRS_struct,nii_name,rot_folder)

%Coregistration of MRS voxel volumes to imaging datasets, based on headers. 

MRS_struct.p.coreg = 1;

if (MRS_struct.ii ~= length(nii_name))
       error('The number of nifti files does not match the number of MRS files processed by GannetLoad.'); 
end

for ii = 1:MRS_struct.ii

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
            if ~exist('rot_folder','var')
                rot_folder=nii_name;
            end
            MRS_struct=GannetMask_GE(fname,nii_name{ii},MRS_struct,rot_folder{ii},ii);
    end
    
    
    
    %Build output figure
        if(ishandle(103))
            close(103)
        end
    
    h=figure(103);
        set(h, 'Position', [100, 100, 1000, 707]);
        set(h,'Color',[1 1 1]);
        figTitle = ['GannetCoRegister Output'];
        set(gcf,'Name',figTitle,'Tag',figTitle, 'NumberTitle','off')

        
         subplot(2,3,4:6)
        axis off;
                tmp = ['Mask output '  ];   
        text(.5,0.75, tmp,'HorizontalAlignment','right',...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        
        
        D = MRS_struct.mask.outfile{ii} ;
        if size(D,2) >30
            [x tmp z] = fileparts(D);
        else
            tmp = [ MRS_struct.mask.outfile(ii) ];
        end
        tmp = regexprep(tmp, '_','-');
        tmp=[': ' tmp]; 
        text(.5,0.75, tmp, ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        
        tmp = ['Spatial Parameters '];
        text(.5,0.63, tmp, 'HorizontalAlignment','right',...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        tmp = [': [LR, AP, FH]'  ];
        text(.5, 0.63, tmp, ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        
        tmp = ['Size ' ];
        text(.5,0.51, tmp, 'HorizontalAlignment','right',...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        tmp = [': [' num2str(MRS_struct.p.voxsize(ii,1)) ', ' num2str(MRS_struct.p.voxsize(ii,2)) ', ' num2str(MRS_struct.p.voxsize(ii,3)) '] mm'  ];
        text(.5,0.51, tmp, ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
 
        
        tmp = ['Volume ' ];
        text(.5,0.39, tmp, 'HorizontalAlignment','right',...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        vol = MRS_struct.p.voxsize(ii,1)*MRS_struct.p.voxsize(ii,2)*MRS_struct.p.voxsize(ii,3)*.001;
        tmp = [': ' num2str(vol) ' ml'  ];
        text(.5, 0.39, tmp, ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
 
        
        tmp = 'Position ';
        text(0.5,0.27, tmp, 'HorizontalAlignment','right',...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        tmp = [': [' num2str(MRS_struct.p.voxoff(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxoff(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxoff(ii,3), '%3.1f') '] mm'  ];
        text(.5, 0.27, tmp, ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);

        
        tmp = 'Angulations ';
        text(0.5,0.15, tmp, 'HorizontalAlignment','right',...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        tmp = [': [' num2str(MRS_struct.p.voxang(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,3), '%3.1f') ', ' '] deg'  ];
        text(.5, 0.15, tmp, ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        
         
        tmp = 'Version (CoRegister) ';
        text(0.5,0.03, tmp, 'HorizontalAlignment','right',...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        tmp = ': 140707';
        text(.5, 0.03, tmp, ...
            'VerticalAlignment', 'top',...
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
        p = get(h,'pos'); % get position of axes
        set(h,'pos',[0.0 0.15 1 1]) % move the axes slightly
        htitle=title(t, 'FontName', 'Helvetica','FontSize',15);
        position = get(htitle,'Position');
        position(2)=position(2)-20;
        set(htitle,'Position',position)       
        script_path=which('GannetLoad');
              Gannet_circle_white=[script_path(1:(end-13)) '/GANNET_circle_white.jpg'];
              A2=imread(Gannet_circle_white);
              hax=axes('Position',[0.80, 0.08, 0.15, 0.15]);
              image(A2);axis off; axis square;

        
       
           %%%%  Save EPS %%%%%
           

epsdirname = [ './MRSCoReg_' datestr(clock,'yymmdd') ];

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
    elseif(strcmpi(MRS_struct.p.vendor,'Philips_data'))  % make this be sdat
        tmp = strfind(pfil_nopath, '.data');
        dot7 = tmp(end); % just in case there's another .data somewhere else...
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
        %epsdirname
        mkdir(epsdirname)
    end
    saveas(gcf, pdfname);
        
        
end