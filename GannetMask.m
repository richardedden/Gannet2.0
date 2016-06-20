function [MRS_struct ] = GannetMask(sparname, nii_file, MRS_struct)



% this relies on SPM, nifti exported by Philips, and spar/sdat

% some code to make ok wth 2 or 3  inputs
% fix nifti inputs and writing output mask to directory of where the nifti
% is. 

[pathspar,namespar,ext] = fileparts(sparname);
[pathnii,namenii,extnii] = fileparts(nii_file);

fidoutmask = fullfile(pathnii,[namespar '_mask.nii'])


%function make_voxel
 
%if nargin =2;   
%     lastchar = sdatfile;   
%     last4char=lastchar((end-3):end);
%     
%     if(strcmp(last4char,'SDAT'))
%        MRS_struct.p.spar_string='SPAR';
%     else
%         MRS_struct.p.spar_string='spar';
%     end
%     
%     sparname = [sdatfile MRS_struct.p.spar_string]
%end
sparheadinfo = textread(sparname, '%s');

sparidx=find(ismember(sparheadinfo, 'ap_size')==1);
MRS_struct.p.voxsize(2) = str2num(sparheadinfo{sparidx+2});
sparidx=find(ismember(sparheadinfo, 'lr_size')==1);
MRS_struct.p.voxsize(1) = str2num(sparheadinfo{sparidx+2});
sparidx=find(ismember(sparheadinfo, 'cc_size')==1);
MRS_struct.p.voxsize(3) = str2num(sparheadinfo{sparidx+2});

sparidx=find(ismember(sparheadinfo, 'ap_off_center')==1);
MRS_struct.p.voxoff(2) = str2num(sparheadinfo{sparidx+2});
sparidx=find(ismember(sparheadinfo, 'lr_off_center')==1);
MRS_struct.p.voxoff(1) = str2num(sparheadinfo{sparidx+2});
sparidx=find(ismember(sparheadinfo, 'cc_off_center')==1);
MRS_struct.p.voxoff(3) = str2num(sparheadinfo{sparidx+2});

sparidx=find(ismember(sparheadinfo, 'ap_angulation')==1);
MRS_struct.p.voxang(2) = str2num(sparheadinfo{sparidx+2});
sparidx=find(ismember(sparheadinfo, 'lr_angulation')==1);
MRS_struct.p.voxang(1) = str2num(sparheadinfo{sparidx+2});
sparidx=find(ismember(sparheadinfo, 'cc_angulation')==1);
MRS_struct.p.voxang(3) = str2num(sparheadinfo{sparidx+2});


V=spm_vol(nii_file);
[T1,XYZ]=spm_read_vols(V);
H=spm_read_hdr(nii_file);

%Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
%tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
halfpixshift = -H.dime.pixdim(1:3).'/2;
halfpixshift(3) = -halfpixshift(3);
XYZ=XYZ+repmat(halfpixshift,[1 size(XYZ,2)]);

% get information from SPAR - change later to be read in

% ap_size = 30;
% lr_size = 30;
% cc_size = 30;
% ap_off = 61.8508873;
% lr_off = -7.173314571;
% cc_off = -0.3336545825;
% ap_ang = -2.264750719;
% lr_ang = -28.66100502;
% cc_ang = 3.604790211;
% 

% ap_size = 40;
% lr_size = 30;
% cc_size = 20;
% ap_off = -3.092470646;
% lr_off = 36.38365555;
% cc_off = -23.80636787;
% ap_ang = 3.970147133;
% lr_ang = 20.98232651;
% cc_ang = 10.02927208;
% 
ap_size = MRS_struct.p.voxsize(2);
lr_size = MRS_struct.p.voxsize(1);
cc_size = MRS_struct.p.voxsize(3);
ap_off = MRS_struct.p.voxoff(2);
lr_off = MRS_struct.p.voxoff(1);
cc_off = MRS_struct.p.voxoff(3);
ap_ang = MRS_struct.p.voxang(2);
lr_ang = MRS_struct.p.voxang(1);
cc_ang = MRS_struct.p.voxang(3);
% 
% 
%We need to flip ap and lr axes to match NIFTI convention
ap_off = -ap_off;
lr_off = -lr_off;


ap_ang = -ap_ang;
lr_ang = -lr_ang;



% define the voxel - use x y z  
% currently have spar convention that have in AUD voxel - will need to
% check for everything in future...
% x - left = positive
% y - posterior = postive
% z - superior = positive
vox_ctr = ...
      [lr_size/2 -ap_size/2 cc_size/2 ;
       -lr_size/2 -ap_size/2 cc_size/2 ;
       -lr_size/2 ap_size/2 cc_size/2 ;
       lr_size/2 ap_size/2 cc_size/2 ;
       -lr_size/2 ap_size/2 -cc_size/2 ;
       lr_size/2 ap_size/2 -cc_size/2 ;
       lr_size/2 -ap_size/2 -cc_size/2 ;
       -lr_size/2 -ap_size/2 -cc_size/2 ];
   
% make rotations on voxel
rad = pi/180;
initrot = zeros(3,3);

xrot = initrot;
xrot(1,1) = 1;
xrot(2,2) = cos(lr_ang *rad);
xrot(2,3) =-sin(lr_ang*rad);
xrot(3,2) = sin(lr_ang*rad);
xrot(3,3) = cos(lr_ang*rad);

yrot = initrot;
yrot(1,1) = cos(ap_ang*rad);
yrot(1,3) = sin(ap_ang*rad);
yrot(2,2) = 1;
yrot(3,1) = -sin(ap_ang*rad);
yrot(3,3) = cos(ap_ang*rad);

zrot = initrot;
zrot(1,1) = cos(cc_ang*rad);
zrot(1,2) = -sin(cc_ang*rad);
zrot(2,1) = sin(cc_ang*rad);
zrot(2,2) = cos(cc_ang*rad);
zrot(3,3) = 1;

% rotate voxel
vox_rot = xrot*yrot*zrot*vox_ctr.';

% calculate corner coordinates relative to xyz origin
vox_ctr_coor = [lr_off ap_off cc_off];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
vox_corner = vox_rot+vox_ctr_coor;



%%% make new comment

mask = zeros(1,size(XYZ,2));
sphere_radius = sqrt((lr_size/2)^2+(ap_size/2)^2+(cc_size/2)^2);
distance2voxctr=sqrt(sum((XYZ-repmat([lr_off ap_off cc_off].',[1 size(XYZ, 2)])).^2,1));
sphere_mask(distance2voxctr<=sphere_radius)=1;
%sphere_mask2=ones(1,(sum(sphere_mask)));

mask(sphere_mask==1) = 1;
XYZ_sphere = XYZ(:,sphere_mask == 1);

tri = delaunayn([vox_corner.'; [lr_off ap_off cc_off]]);
tn = tsearchn([vox_corner.'; [lr_off ap_off cc_off]], tri, XYZ_sphere.');
isinside = ~isnan(tn);
mask(sphere_mask==1) = isinside;

mask = reshape(mask, V.dim);

V_mask.fname=[fidoutmask];
V_mask.descrip='MRS_Voxel_Mask';
V_mask.dim=V.dim;
V_mask.dt=V.dt;
V_mask.mat=V.mat;

V_mask=spm_write_vol(V_mask,mask);

T1img = T1/max(T1(:));
T1img_mas = T1img + .2*mask;

% construct output
% 
 voxel_ctr = [-lr_off -ap_off cc_off];



%MRS_struct.mask.dim(MRS_struct.ii,:)=V.dim;
%MRS_struct.mask.img(MRS_struct.ii,:,:,:)=T1img_mas;
MRS_struct.mask.outfile(MRS_struct.ii,:)=fidoutmask;
%This assumes 1-mm iso T1 - need to fix at a later date.
slice = [round(V.dim(1)/2+voxel_ctr(1)) 
        round(V.dim(2)/2+voxel_ctr(2)) 
        round(V.dim(3)/2+voxel_ctr(3))];

size_max=max(size(T1img_mas));
three_plane_img=zeros([size_max 3*size_max]);
im1 = squeeze(T1img_mas(:,:,slice(3)));
im1 = im1(end:-1:1,:)';  %not sure if need this '
im3 = squeeze(T1img_mas(:,slice(2),:));
im3 = im3(end:-1:1,end:-1:1)'; %may not need '
im2 = squeeze(T1img_mas(slice(1),:,:));
im2 = im2(:,end:-1:1)';

three_plane_img(:,1:size_max) = image_center(im1, size_max);
three_plane_img(:,size_max*2+(1:size_max))=image_center(im3,size_max);
three_plane_img(:,size_max+(1:size_max))=image_center(im2,size_max);

MRS_struct.mask.img(MRS_struct.ii,:,:)=three_plane_img;


%figure(198)
%imagesc(three_plane_img);
%colormap('gray');
%caxis([0 1])
%axis equal;
%axis tight;
%axis off;

%end

