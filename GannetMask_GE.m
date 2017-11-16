function MRS_struct = GannetMask_GE(Pname, dcm_dir, MRS_struct, dcm_dir2,ii)

if nargin == 3
    if isstruct(MRS_struct)
        %dcm_dir2={dcm_dir};
        dcm_dir2=dcm_dir; %ADH - don't think should have {} but might break Gannet batching for GABA
        ii =1;
    else
        dcm_dir2=MRS_struct;
        clear MRS_struct;
        MRS_struct.ii=1;
        ii = 1;
    end
    
end

if nargin == 4
    if isnumeric(dcm_dir2)
        ii = dcm_dir2;
        dcm_dir2 = dcm_dir;
    else
        ii = 1;
    end
end

if nargin < 3
    MRS_struct.ii=1;
    ii = 1;
    %dcm_dir2={dcm_dir};
    dcm_dir2=dcm_dir; %ADH - don't think should have {} but might break Gannet batching for GABA
end

% This relies on SPM and anatomical images as dicoms

% Parse P-file to extract voxel dimensions and offsets (MM: 171110)
if MRS_struct.p.GE.rdbm_rev_num >= 11.0
    % In 11.0 and later the header and data are stored as little-endian
    fid = fopen(Pname,'r', 'ieee-le');
else
    fid = fopen(Pname, 'r', 'ieee-be');
end
fseek(fid, 1468, 'bof');
p_hdr_value = fread(fid, 12, 'integer*4'); % byte offsets to start of sub-header structures
fseek(fid, p_hdr_value(8), 'bof'); % set position to start of rdb_hdr_exam
o_hdr_value = fread(fid, p_hdr_value(9)-p_hdr_value(8), 'real*4');
if MRS_struct.p.GE.rdbm_rev_num == 16
    MRS_struct.p.voxdim(ii,:) = o_hdr_value(822:824)';
    MRS_struct.p.voxoff(ii,:) = o_hdr_value(825:827)';
elseif MRS_struct.p.GE.rdbm_rev_num == 24
    MRS_struct.p.voxdim(ii,:) = o_hdr_value(1228:1230)';
    MRS_struct.p.voxoff(ii,:) = o_hdr_value(1231:1233)';
end
MRS_struct.p.voxoff(ii,:) = MRS_struct.p.voxoff(ii,:) .* [-1 -1 1];

% ptr = fopen([Pname '.hdr']);
% MRSHead = fscanf(ptr,'%c');
% fclose(ptr);

% k = findstr(' user11:', MRSHead);
% b = str2double(MRSHead(k+11:k+18));
% k = findstr(' user12:', MRSHead);
% c = str2double(MRSHead(k+11:k+18));
% k = findstr(' user13:', MRSHead);
% a = str2double(MRSHead(k+11:k+18));
% MRS_struct.p.voxoff(ii,:) = [-b -c a];

% k = findstr(' user8:', MRSHead);
% d = str2double(MRSHead(k+10:k+17));
% k = findstr(' user9:', MRSHead);
% e = str2double(MRSHead(k+10:k+17));
% %e = 0.5*e; %e is the anterior posterior direction
% k = findstr(' user10:', MRSHead);
% f = str2double(MRSHead(k+11:k+18));
% MRS_struct.p.voxsize(ii,:) = [d e f]; % works for ob-axial rotator

% MRS_struct.p.voxang is not contained in P-file header (really!)
% The rotation is adopted from the image on which the voxel was placed
% i.e. either the 3D T1 or a custom rotated localizer.
MRS_struct.p.voxang(ii,:) = [NaN NaN NaN]; % put as NaN for now - for output page
currdir = pwd;

cd(dcm_dir2);
dcm_list = dir;

for k = length(dcm_list):-1:1
    fname = dcm_list(k).name;
    if fname(1) == '.'
        dcm_list(k) = [];
    end
end

dcm_list2 = dcm_list.name;
if dcm_list2((end-2:end))=='nii'
    dcm_list2 = dcm_list(4).name;
end

% Load dicoms
MRSRotHead = dicominfo(dcm_list2);

cd(currdir);
%Do some housekeeping on dicom series
currdir=pwd;
cd(dcm_dir);
dcm_list=dir;
dcm_list=dcm_list(3:end);
hdr = spm_dicom_headers(char(dcm_list.name));
spm_dicom_convert(hdr, 'all','flat','nii');
nii_file_dir=dir('s*.nii');
cd(currdir);

nii_file=[dcm_dir '/' nii_file_dir(1).name];

V=spm_vol(nii_file);
[T1,XYZ]=spm_read_vols(V);
H=spm_read_hdr(nii_file);

%Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
%tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
halfpixshift = -H.dime.pixdim(1:3).'/2;
%halfpixshift(3) = -halfpixshift(3);
XYZ=XYZ+repmat(halfpixshift,[1 size(XYZ,2)]);

ap_size = MRS_struct.p.voxdim(ii,2);
lr_size = MRS_struct.p.voxdim(ii,1);
cc_size = MRS_struct.p.voxdim(ii,3);
ap_off = MRS_struct.p.voxoff(ii,2);
lr_off = MRS_struct.p.voxoff(ii,1);
cc_off = MRS_struct.p.voxoff(ii,3);
%ap_ang = MRS_struct.p.voxang(2);
%lr_ang = MRS_struct.p.voxang(1);
%cc_ang = MRS_struct.p.voxang(3);

%We need to flip ap and lr axes to match NIFTI convention
ap_off = -ap_off;
lr_off = -lr_off;

% define the voxel - use x y z
% x - left = positive
% y - posterior = postive
% z - superior = positive
% vox_ctr = ...
%     [lr_size/2 -ap_size/2 cc_size/2 ;
%     -lr_size/2 -ap_size/2 cc_size/2 ;
%     -lr_size/2 ap_size/2 cc_size/2 ;
%     lr_size/2 ap_size/2 cc_size/2 ;
%     -lr_size/2 ap_size/2 -cc_size/2 ;
%     lr_size/2 ap_size/2 -cc_size/2 ;
%     lr_size/2 -ap_size/2 -cc_size/2 ;
%     -lr_size/2 -ap_size/2 -cc_size/2 ];

%vox_rot=rotmat*vox_ctr.';

% calculate corner coordinates relative to xyz origin
vox_ctr_coor = [lr_off ap_off cc_off];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
% vox_corner = vox_rot+vox_ctr_coor;

%New code RAEE

MRS_Rot_RE=MRSRotHead.ImageOrientationPatient;
MRS_Rot_RE=reshape(MRS_Rot_RE',[3 2]);
edge1=repmat(MRS_Rot_RE(:,1), [1 8]);
edge1(:,2:3)=-edge1(:,2:3);
edge1(:,5)=-edge1(:,5);
edge1(:,8)=-edge1(:,8);
edge1(1:2,:)=-edge1(1:2,:);

edge2=repmat(MRS_Rot_RE(:,2), [1 8]);
edge2(:,1:2)=-edge2(:,1:2);
edge2(:,7)=-edge2(:,7);
edge2(:,8)=-edge2(:,8);
edge2(1:2,:)=-edge2(1:2,:);
edge3=repmat(cross(MRS_Rot_RE(:,1),MRS_Rot_RE(:,2)),[1 8]);
edge3(:,5:8)=-edge3(:,5:8);
edge3(1:2,:)=-edge3(1:2,:);
vox_corner=vox_ctr_coor+lr_size/2*edge1+ap_size/2*edge2+cc_size/2*edge3;

mask = zeros(1,size(XYZ,2));
sphere_radius = sqrt((lr_size/2)^2+(ap_size/2)^2+(cc_size/2)^2);
distance2voxctr = sqrt(sum((XYZ-repmat([lr_off ap_off cc_off].',[1 size(XYZ, 2)])).^2,1));
sphere_mask(distance2voxctr <= sphere_radius) = 1;

mask(sphere_mask==1) = 1;
XYZ_sphere = XYZ(:,sphere_mask == 1);

tri = delaunayn([vox_corner.'; [lr_off ap_off cc_off]]);
tn = tsearchn([vox_corner.'; [lr_off ap_off cc_off]], tri, XYZ_sphere.');
isinside = ~isnan(tn);
mask(sphere_mask==1) = isinside;

mask = reshape(mask, V.dim);

[~,namespar] = fileparts(Pname);
fidoutmask = fullfile(dcm_dir,[namespar '_mask.nii']);
V_mask.fname=fidoutmask;
V_mask.descrip='MRS_Voxel_Mask';
V_mask.dim=V.dim;
V_mask.dt=V.dt;
V_mask.mat=V.mat;

spm_write_vol(V_mask,mask);

T1img = T1/max(T1(:));
T1img_mas = T1img + .2*mask;

% construct output

fidoutmask = cellstr(fidoutmask);
MRS_struct.mask.outfile(ii,:) = fidoutmask;

voxel_ctr = [-lr_off -ap_off cc_off];
voxel_ctr(1:2)=-voxel_ctr(1:2);
voxel_search=(XYZ(:,:)-repmat(voxel_ctr.',[1 size(XYZ,2)])).^2;
voxel_search=sqrt(sum(voxel_search,1));
[~,index1]=min(voxel_search);
[slice(1), slice(2), slice(3)]=ind2sub( V.dim,index1);

size_max=max(size(T1img_mas));
three_plane_img=zeros([size_max 3*size_max]);
im1 = squeeze(T1img_mas(:,:,slice(3)));
im1 = im1(end:-1:1,end:-1:1)';  %not sure if need this '
im3 = squeeze(T1img_mas(:,slice(2),:));
im3 = im3(end:-1:1,end:-1:1)'; %may not need '
im2 = squeeze(T1img_mas(slice(1),:,:));
im2 = im2(end:-1:1,end:-1:1)';

three_plane_img(:,1:size_max) = image_center(im1, size_max);
three_plane_img(:,size_max*2+(1:size_max))=image_center(im3,size_max);
three_plane_img(:,size_max+(1:size_max))=image_center(im2,size_max);

MRS_struct.mask.img(ii,:,:)=three_plane_img;
MRS_struct.mask.T1image(ii,:) = {nii_file};

% figure(198)
% imagesc(three_plane_img);
% colormap('gray');
% caxis([0 1]);
% axis equal;
% axis tight;
% axis off;

end

