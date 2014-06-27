
function MRS_struct = GannetCoRegister(MRS_struct,nii_name)

%Coregistration of MRS voxel volumes to imaging datasets, based on headers. 

MRS_struct.p.coreg = 1;
%Ultimately this switch will not be necessary...
    switch MRS_struct.p.vendor
    
    case 'Philips'

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
        error(['GannetCoRegister does not yet support ' MRS_struct.p.vendor ' data.']);
    end
    
    if (MRS_struct.ii ~= length(nii_name))
       error('The number of nifti files does not match the number of MRS files processed by GannetLoad.'); 
    end
    %Currently only SDAT is supported
    %Run the script...
    for ii=1:length(nii_name)
        fname = MRS_struct.gabafile{ii};
        sparname = [fname(1:(end-4)) MRS_struct.p.spar_string];
        [MRS_struct Mask]=GannetMask(sparname,nii_name{ii},MRS_struct);
    end
    
    %Build output figure
    h=figure(103);
        set(h, 'Position', [100, 100, 1000, 707]);
        set(h,'Color',[1 1 1]);
        figTitle = ['GannetCoRegister Output'];
        set(gcf,'Name',figTitle,'Tag',figTitle, 'NumberTitle','off');
              
        
        voxel_ctr = MRS_struct.p.voxoff;
        voxel_ctr(1:2)=-voxel_ctr(1:2);
        
        
        %This assumes 1-mm iso T1 - need to fix at a later date.
        slice = [round(Mask.dim(1)/2+voxel_ctr(1)) 
                round(Mask.dim(2)/2+voxel_ctr(2)) 
                round(Mask.dim(3)/2+voxel_ctr(3))];

        size_max=max(size(Mask.img));
        three_plane_img=zeros([size_max 3*size_max]);
        im1 = squeeze(Mask.img(:,:,slice(3)));
        im1 = im1';  %not sure if need this '
        im2 = squeeze(Mask.img(:,slice(2),:));
        im2 = im2(:,end:-1:1)'; %may not need '
        im3 = squeeze(Mask.img(slice(1),:,:));
        im3 = im3(:,end:-1:1)';

        three_plane_img(:,1:size_max) = image_center(im1, size_max);
        three_plane_img(:,size_max*2+(1:size_max))=image_center(im2,size_max);
        three_plane_img(:,size_max+(1:size_max))=image_center(im3,size_max);
        h=subplot(2,2,1:2);
        imagesc(three_plane_img);
        colormap('gray');
        caxis([0 1])
        axis equal;
        axis tight;
        axis off;
        p = get(h,'pos') % get position of axes
        set(h,'pos',[0.1 0.5 0.8 0.4]) % move the axes slightly
        
        script_path=which('GannetLoad');
              % CJE update for GE
    %          Gannet_circle=[script_path(1:(end-12)) 'GANNET_circle.png'];
              Gannet_circle_white=[script_path(1:(end-13)) '/GANNET_circle_white.jpg'];
    %          A=imread(Gannet_circle);
              A2=imread(Gannet_circle_white);
              hax=axes('Position',[0.80, 0.05, 0.15, 0.15]);
              %set(gca,'Units','normalized');set(gca,'Position',[0.05 0.05 1.85 0.15]);
              image(A2);axis off; axis square;

        %axis off;

end