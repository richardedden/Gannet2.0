function MRS_struct = GannetMask_Philips(sparname, nii_file, MRS_struct, ii, vox, kk)

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:qhullmx:InternalWarning');

if nargin == 2
    MRS_struct.ii = 1;
    ii = 1;
end

[~,metabfile] = fileparts(sparname);
pathnii = fileparts(nii_file);

fidoutmask = fullfile(pathnii,[metabfile '_mask.nii']);

sparname = fopen(sparname,'r');
sparheadinfo = textscan(sparname, '%s');
sparheadinfo = sparheadinfo{1};

sparidx = find(ismember(sparheadinfo, 'ap_size')==1);
MRS_struct.p.voxdim(ii,2) = str2double(sparheadinfo{sparidx+2});
sparidx = find(ismember(sparheadinfo, 'lr_size')==1);
MRS_struct.p.voxdim(ii,1) = str2double(sparheadinfo{sparidx+2});
sparidx = find(ismember(sparheadinfo, 'cc_size')==1);
MRS_struct.p.voxdim(ii,3) = str2double(sparheadinfo{sparidx+2});

sparidx = find(ismember(sparheadinfo, 'ap_off_center')==1);
MRS_struct.p.voxoff(ii,2) = str2double(sparheadinfo{sparidx+2});
sparidx = find(ismember(sparheadinfo, 'lr_off_center')==1);
MRS_struct.p.voxoff(ii,1) = str2double(sparheadinfo{sparidx+2});
sparidx = find(ismember(sparheadinfo, 'cc_off_center')==1);
MRS_struct.p.voxoff(ii,3) = str2double(sparheadinfo{sparidx+2});

sparidx = find(ismember(sparheadinfo, 'ap_angulation')==1);
MRS_struct.p.voxang(ii,2) = str2double(sparheadinfo{sparidx+2});
sparidx = find(ismember(sparheadinfo, 'lr_angulation')==1);
MRS_struct.p.voxang(ii,1) = str2double(sparheadinfo{sparidx+2});
sparidx = find(ismember(sparheadinfo, 'cc_angulation')==1);
MRS_struct.p.voxang(ii,3) = str2double(sparheadinfo{sparidx+2});

V = spm_vol(nii_file);
[T1,XYZ] = spm_read_vols(V);

% Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
% tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
[~,voxdim] = spm_get_bbox(V,'fv'); % MM (180220)
voxdim = abs(voxdim)';
halfpixshift = -voxdim(1:3)/2;
halfpixshift(3) = -halfpixshift(3);
XYZ = XYZ + repmat(halfpixshift, [1 size(XYZ,2)]);

% Get information from SPAR - change later to be read in
ap_size = MRS_struct.p.voxdim(ii,2);
lr_size = MRS_struct.p.voxdim(ii,1);
cc_size = MRS_struct.p.voxdim(ii,3);
ap_off  = MRS_struct.p.voxoff(ii,2);
lr_off  = MRS_struct.p.voxoff(ii,1);
cc_off  = MRS_struct.p.voxoff(ii,3);
ap_ang  = MRS_struct.p.voxang(ii,2);
lr_ang  = MRS_struct.p.voxang(ii,1);
cc_ang  = MRS_struct.p.voxang(ii,3);

% We need to flip ap and lr axes to match NIFTI convention
ap_off = -ap_off;
lr_off = -lr_off;
ap_ang = -ap_ang;
lr_ang = -lr_ang;

% Define the voxel - use x y z
% Currently have spar convention that have in AUD voxel - will need to
% check for everything in future...
% x - left = positive
% y - posterior = postive
% z - superior = positive
vox_ctr = ...
    [lr_size/2 -ap_size/2  cc_size/2;
    -lr_size/2 -ap_size/2  cc_size/2;
    -lr_size/2  ap_size/2  cc_size/2;
     lr_size/2  ap_size/2  cc_size/2;
    -lr_size/2  ap_size/2 -cc_size/2;
     lr_size/2  ap_size/2 -cc_size/2;
     lr_size/2 -ap_size/2 -cc_size/2;
    -lr_size/2 -ap_size/2 -cc_size/2];

% Make rotations on voxel
rad = pi/180;
initrot = zeros(3,3);

xrot      = initrot;
xrot(1,1) = 1;
xrot(2,2) = cos(lr_ang *rad);
xrot(2,3) = -sin(lr_ang*rad);
xrot(3,2) = sin(lr_ang*rad);
xrot(3,3) = cos(lr_ang*rad);

yrot      = initrot;
yrot(1,1) = cos(ap_ang*rad);
yrot(1,3) = sin(ap_ang*rad);
yrot(2,2) = 1;
yrot(3,1) = -sin(ap_ang*rad);
yrot(3,3) = cos(ap_ang*rad);

zrot      = initrot;
zrot(1,1) = cos(cc_ang*rad);
zrot(1,2) = -sin(cc_ang*rad);
zrot(2,1) = sin(cc_ang*rad);
zrot(2,2) = cos(cc_ang*rad);
zrot(3,3) = 1;

% Rotate voxel
vox_rot = xrot * yrot * zrot * vox_ctr.';

% Calculate corner coordinates relative to xyz origin
vox_ctr_coor = [lr_off ap_off cc_off];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
vox_corner = vox_rot+vox_ctr_coor;

mask = zeros(1,size(XYZ,2));
sphere_radius = sqrt((lr_size/2)^2+(ap_size/2)^2+(cc_size/2)^2);
distance2voxctr = sqrt(sum((XYZ-repmat([lr_off ap_off cc_off].',[1 size(XYZ,2)])).^2,1));
sphere_mask(distance2voxctr <= sphere_radius) = 1;

mask(sphere_mask == 1) = 1;
XYZ_sphere = XYZ(:,sphere_mask == 1);

tri = delaunayn([vox_corner.'; [lr_off ap_off cc_off]]);
tn = tsearchn([vox_corner.'; [lr_off ap_off cc_off]], tri, XYZ_sphere.');
isinside = ~isnan(tn);
mask(sphere_mask==1) = isinside;

mask = reshape(mask, V.dim);

V_mask.fname   = fidoutmask;
V_mask.descrip = 'MRS_voxel_mask';
V_mask.dim     = V.dim;
V_mask.dt      = V.dt;
V_mask.mat     = V.mat;

V_mask = spm_write_vol(V_mask,mask);

% Build output

fidoutmask = cellstr(fidoutmask);
MRS_struct.mask.(vox{kk}).outfile(ii,:) = fidoutmask;

voxel_ctr      = [-lr_off -ap_off cc_off];
voxel_ctr(1:2) = -voxel_ctr(1:2);

% Transform structural image and co-registered voxel mask from voxel to
% world space for output (MM: 180221)
[img_t,img_c,img_s] = voxel2world_space(V,voxel_ctr);
[mask_t,mask_c,mask_s] = voxel2world_space(V_mask,voxel_ctr);

img_t = flipud(img_t/max(T1(:)));
img_c = flipud(img_c/max(T1(:)));
img_s = flipud(img_s/max(T1(:)));

img_t = img_t + 0.2*flipud(mask_t);
img_c = img_c + 0.2*flipud(mask_c);
img_s = img_s + 0.2*flipud(mask_s);

size_max = max([max(size(img_t)) max(size(img_c)) max(size(img_s))]);
three_plane_img = zeros([size_max 3*size_max]);
three_plane_img(:,1:size_max)              = image_center(img_t, size_max);
three_plane_img(:,size_max+(1:size_max))   = image_center(img_s, size_max);
three_plane_img(:,size_max*2+(1:size_max)) = image_center(img_c, size_max);

MRS_struct.mask.(vox{kk}).img{ii} = three_plane_img;
MRS_struct.mask.(vox{kk}).T1image(ii,:) = {nii_file};

warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:qhullmx:InternalWarning');

end



