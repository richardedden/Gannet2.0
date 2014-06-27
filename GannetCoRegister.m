
function MRS_struct = GannetCoRegister(MRS_struct,nii_name)

%Coregistration of MRS voxel volumes to imaging datasets, based on headers. 

MRS_struct.p.coreg = 1;

if (MRS_struct.ii ~= length(nii_name))
       error('The number of nifti files does not match the number of MRS files processed by GannetLoad.'); 
end

%Ultimately this switch will not be necessary...
    switch MRS_struct.p.vendor
    
    case 'Philips'
        for ii=1:length(nii_name)
            fname = MRS_struct.gabafile{ii};
            sparname = [fname(1:(end-4)) MRS_struct.p.spar_string];
            MRS_struct=GannetMask_Philips(sparname,nii_name{ii},MRS_struct);
        end
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
        for ii=1:length(nii_name)
            fname = MRS_struct.gabafile{ii};
            %sparname = [fname(1:(end-4)) MRS_struct.p.spar_string];
            MRS_struct=GannetMask_GE(fname,nii_name{ii},MRS_struct);
        end
        %error(['GannetCoRegister does not yet support ' MRS_struct.p.vendor ' data.']);
    end
    
    
    %Currently only SDAT is supported
    %Run the script...
%    for ii=1:length(nii_name)
%        fname = MRS_struct.gabafile{ii};
%        sparname = [fname(1:(end-4)) MRS_struct.p.spar_string];
%        MRS_struct=GannetMask(sparname,nii_name{ii},MRS_struct);
%    end
    
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

end