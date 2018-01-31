function MRS_struct = GannetMask_SiemensTWIX(filename, nii_file, MRS_struct, ii)

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:qhullmx:InternalWarning');

% being testing on data from Univ of Florida - will need to extend
% also may need to extend to change dicoms into nifti

[path,name] = fileparts(filename);
fidoutmask = fullfile(path,[name '_mask.nii']);

fid = fopen(filename);
line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_Normal_Sag">  { <Precision> ');
equals_index = strfind(line,'16 ');
while isempty(index) || isempty(equals_index)
    line = fgets(fid);
    index = strfind(line,'<ParamDouble."VoI_Normal_Sag">  { <Precision> ');
    equals_index = strfind(line,'16 ');
end
% GO 07/14/16: If a parameter is set to zero (e.g. if no voxel rotation is
% performed), the respective field is left empty in the TWIX file. This
% case needs to be intercepted. Same goes for all parameters below that are
% extracted from the TWIX file. Setting it to 0 may throw up calculations
% below, so set to the minimum possible.
if strcmp(line(equals_index+4),'}')
    NormSag = 0;
else
    NormSag = line(equals_index+4:equals_index+16);
    NormSag = str2double(NormSag);
end

line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_Normal_Cor">  { <Precision> ');
equals_index = strfind(line,'16 ');
while isempty(index) || isempty(equals_index)
    line = fgets(fid);
    index = strfind(line,'<ParamDouble."VoI_Normal_Cor">  { <Precision> ');
    equals_index = strfind(line,'16 ');
end
if strcmp(line(equals_index+4),'}')
    NormCor = 0;
else
    NormCor = line(equals_index+4:equals_index+16);
    NormCor = str2double(NormCor);
end

line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_Normal_Tra">  { <Precision> ');
equals_index = strfind(line,'16 ');
while isempty(index) || isempty(equals_index)
    line = fgets(fid);
    index = strfind(line,'<ParamDouble."VoI_Normal_Tra">  { <Precision> ');
    equals_index = strfind(line,'16 ');
end
if strcmp(line(equals_index+4),'}')
    NormTra = 0;
else
    NormTra = line(equals_index+4:equals_index+16);
    NormTra = str2double(NormTra);
end

line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_Position_Sag">  { <Precision> ');
equals_index = strfind(line,'16 ');
while isempty(index) || isempty(equals_index)
    line=fgets(fid);
    index = strfind(line,'<ParamDouble."VoI_Position_Sag">  { <Precision> ');
    equals_index = strfind(line,'16 ');
end
if strcmp(line(equals_index+4),'}')
    PosSag = realmin('double');
else
    PosSag = line(equals_index+4:equals_index+16);
    PosSag = str2double(PosSag);
end

line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_Position_Cor">  { <Precision> ');
equals_index = strfind(line,'16 ');
while isempty(index) || isempty(equals_index)
    line = fgets(fid);
    index = strfind(line,'<ParamDouble."VoI_Position_Cor">  { <Precision> ');
    equals_index = strfind(line,'16 ');
end
if strcmp(line(equals_index+4),'}')
    PosCor = realmin('double');
else
    PosCor = line(equals_index+4:equals_index+16);
    PosCor = str2double(PosCor);
end

line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_Position_Tra">  { <Precision> ');
equals_index = strfind(line,'16 ');
while isempty(index) || isempty(equals_index)
    line = fgets(fid);
    index = strfind(line,'<ParamDouble."VoI_Position_Tra">  { <Precision> ');
    equals_index = strfind(line,'16 ');
end
if strcmp(line(equals_index+4),'}')
    PosTra = realmin('double');
else
    PosTra = line(equals_index+4:equals_index+16);
    PosTra = str2double(PosTra);
end

line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_SliceThickness">  { <Precision> ');
equals_index=strfind(line,'16 ');
while isempty(index) || isempty(equals_index)
    line=fgets(fid);
    index=strfind(line,'<ParamDouble."VoI_SliceThickness">  { <Precision> ');
    equals_index=strfind(line,'16 ');
end
if strcmp(line(equals_index+4),'}')
    VOIThickness = realmin('double');
else
    VOIThickness = line(equals_index+4:equals_index+16);
    VOIThickness = str2double(VOIThickness);
end

fclose(fid);
fid = fopen(filename);
line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_InPlaneRotAngle">  { <Precision> ');
equals_index = strfind(line,'16 ');
while isempty(index) || isempty(equals_index)
    line = fgets(fid);
    index = strfind(line,'<ParamDouble."VoI_InPlaneRotAngle">  { <Precision> ');
    equals_index = strfind(line,'16 ');
end
if strcmp(line(equals_index+4),'}')
    VoI_InPlaneRot = realmin('double');
else
    VoI_InPlaneRot = line(equals_index+4:equals_index+16);
    VoI_InPlaneRot = str2double(VoI_InPlaneRot);
end

line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_RoFOV">  { <Precision> ');
equals_index = strfind(line,'16 ');
while isempty(index) || isempty(equals_index)
    line = fgets(fid);
    index = strfind(line,'<ParamDouble."VoI_RoFOV">  { <Precision> ');
    equals_index = strfind(line,'16 ');
end
if strcmp(line(equals_index+4),'}')
    VoI_RoFOV = realmin('double');
else
    VoI_RoFOV = line(equals_index+4:equals_index+16);
    VoI_RoFOV = str2double(VoI_RoFOV);
end

line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_PeFOV">  { <Precision> ');
equals_index=strfind(line,'16 ');
while isempty(index) || isempty(equals_index)
    line = fgets(fid);
    index = strfind(line,'<ParamDouble."VoI_PeFOV">  { <Precision> ');
    equals_index = strfind(line,'16 ');
end
if strcmp(line(equals_index+4),'}')
    VoI_PeFOV = realmin('double');
else
    VoI_PeFOV = line(equals_index+4:equals_index+16);
    VoI_PeFOV = str2double(VoI_PeFOV);
end

fclose(fid);

% Beging voxel co-registration

ZED = [NormSag NormCor NormTra];
ZED = ZED *-1;
ROT = VoI_InPlaneRot;

R(1,1) = cos(ROT)+(ZED(1)^2) * (1-cos(ROT));
R(1,2) = ZED(1)*ZED(2)*(1-cos(ROT))-(ZED(3)*sin(ROT));
R(1,3) = ZED(1)*ZED(3)*(1-cos(ROT))+(ZED(2)*sin(ROT));
R(2,1) = ZED(2)*ZED(1)*(1-cos(ROT))+(ZED(3)*sin(ROT));
R(2,2) = cos(ROT)+(ZED(2)^2) * (1-cos(ROT));
R(2,3) = ZED(2)*ZED(3)*(1-cos(ROT))-(ZED(1)*sin(ROT));
R(3,1) = ZED(3)*ZED(1)*(1-cos(ROT))-(ZED(2)*sin(ROT));
R(3,2) = ZED(3)*ZED(2)*(1-cos(ROT))+(ZED(1)*sin(ROT));
R(3,3) = cos(ROT)+(ZED(3)^2) * (1-cos(ROT));

Raxisangle = R;
Raxisangle(1,4) = 0;
Raxisangle(2,4) = 0;
Raxisangle(3,4) = 0;
Raxisangle(4,:) = [0 0 0 1];

colStart(1,1) = 0;
colStart(3,1) = 1/sqrt(((-1.0 * ZED(3) / ZED(2))^2)+1);
colStart(2,1) = -1.0 * colStart(3) * ZED(3) / ZED(2);

ColStart = colStart;
ColStart(1:3,4) = 0;
ColStart(4,:) = [0 0 0 1];

ColMat = Raxisangle*ColStart;
Col(1) = ColMat(1,1);
Col(2) = ColMat(2,1);
Col(3) = ColMat(3,1);

%Col

Row(1) = Col(2) * ZED(3) - Col(3) * ZED(2);
Row(2) = Col(3) * ZED(1) - Col(1) * ZED(3);
Row(3) = Col(1) * ZED(2) - Col(2) * ZED(1);

%Row

%MRS_struct.p.voxoff=[ rda.position(1) rda.position(2) rda.position(3)];
MRS_struct.p.voxoff(ii,:) = [PosSag PosCor PosTra];
MRS_struct.p.voxdim(ii,:) = [VoI_RoFOV VoI_PeFOV VOIThickness];
MRS_Rot(:,1) = Row';
MRS_Rot(:,2) = Col';

MRS_Rot(:,1) = MRS_Rot(:,1) .* [-1 -1 1]';
MRS_Rot(:,2) = MRS_Rot(:,2) .* [-1 -1 1]';

MRS_Rot(:,3) = cross(MRS_Rot(:,2),MRS_Rot(:,1));

rotmat=MRS_Rot; % used to be rotmat = -MRS_rot but doesn't seem to do anything
%- bc of corners defn in mat - the values just get reorded in the vox_corner mat
%rotmat(:,3) = [1 1 1];
%error('gothere');

% voxang is initialized to zero so the code runs, but ideally, it needs to
% parse the rotation info from nii_file2.

V = spm_vol(nii_file);
[T1,XYZ] = spm_read_vols(V);

%Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
%tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
voxdim = abs(V.mat(1:3,1:3));
voxdim = voxdim(voxdim > 0.1);
voxdim = round(voxdim/1e-2)*1e-2; % MM (180130)
halfpixshift = -voxdim(1:3)/2;
halfpixshift(3) = -halfpixshift(3);
XYZ = XYZ + repmat(halfpixshift, [1 size(XYZ,2)]);

ap_size = MRS_struct.p.voxdim(ii,2);
lr_size = MRS_struct.p.voxdim(ii,1);
cc_size = MRS_struct.p.voxdim(ii,3);
ap_off = MRS_struct.p.voxoff(ii,2);
lr_off = MRS_struct.p.voxoff(ii,1);
cc_off = MRS_struct.p.voxoff(ii,3);
%ap_ang = MRS_struct.p.voxang(2);
%lr_ang = MRS_struct.p.voxang(1);
%cc_ang = MRS_struct.p.voxang(3);

% We need to flip ap and lr axes to match NIFTI convention
ap_off = -ap_off;
lr_off = -lr_off;
%ap_ang = -ap_ang;
%lr_ang = -lr_ang;

% define the voxel - use xyz
% currently have spar convention - will need to
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

% vox_ctr = ...
%       [lr_size 0 cc_size ;
%       0 0 cc_size ;
%        0 ap_size cc_size ;
%        lr_size ap_size cc_size ;
%        0 ap_size 0 ;
%        lr_size ap_size 0 ;
%        lr_size 0 0 ;
%        0 0 0];

vox_rot = rotmat*vox_ctr.';

% Calculate corner coordinates relative to xyz origin
vox_ctr_coor = [lr_off ap_off cc_off];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
vox_corner = vox_rot + vox_ctr_coor;

mask = zeros(1,size(XYZ,2));
sphere_radius = sqrt((lr_size/2)^2+(ap_size/2)^2+(cc_size/2)^2);
distance2voxctr = sqrt(sum((XYZ-repmat([lr_off ap_off cc_off].',[1 size(XYZ, 2)])).^2,1));
sphere_mask(distance2voxctr <= sphere_radius) = 1;

mask(sphere_mask == 1) = 1;
XYZ_sphere = XYZ(:,sphere_mask == 1);

tri = delaunayn([vox_corner.'; [lr_off ap_off cc_off]]);
tn = tsearchn([vox_corner.'; [lr_off ap_off cc_off]], tri, XYZ_sphere.');
isinside = ~isnan(tn);
mask(sphere_mask==1) = isinside;

mask = reshape(mask, V.dim);

V_mask.fname   = fidoutmask ;
V_mask.descrip = 'MRS_voxel_mask';
V_mask.dim     = V.dim;
V_mask.dt      = V.dt;
V_mask.mat     = V.mat;

spm_write_vol(V_mask,mask);

T1img = T1/max(T1(:));
T1img_mas = T1img + 0.2*mask;

% Build output

%FOR NOW NEED TO FIX
fidoutmask = cellstr(fidoutmask);
MRS_struct.mask.outfile(ii,:) = fidoutmask;
MRS_struct.p.voxang(ii,:) = [NaN NaN NaN];  % put as NaN for now - for output page
% this is similar to GE - don't have the angles directly - can fix later

voxel_ctr      = [-lr_off -ap_off cc_off];
voxel_ctr(1:2) = -voxel_ctr(1:2);
voxel_search   = (XYZ(:,:)-repmat(voxel_ctr.',[1 size(XYZ,2)])).^2;
voxel_search   = sqrt(sum(voxel_search,1));
[~,index1]     = min(voxel_search);
[slice(1), slice(2), slice(3)] = ind2sub(V.dim,index1);

im1 = squeeze(T1img_mas(:,:,slice(3)));
im1 = im1(end:-1:1,end:-1:1)';
im3 = squeeze(T1img_mas(:,slice(2),:));
im3 = im3(end:-1:1,end:-1:1)';
im2 = squeeze(T1img_mas(slice(1),:,:));
im2 = im2(end:-1:1,end:-1:1)';

% MM (180130): Resize slices if voxel resolution in T1 image isn't
% isometric
if voxdim(1) ~= voxdim(2)
    a = max(voxdim([1 2])) ./ min(voxdim([1 2]));
    im1 = imresize(im1, [size(im1,1)*a size(im1,2)]);
end

if voxdim(1) ~= voxdim(3)
    a = max(voxdim([1 3])) ./ min(voxdim([1 3]));
    im3 = imresize(im3, [size(im3,1)*a size(im3,2)]);
end

if voxdim(2) ~= voxdim(3)
    a = max(voxdim([2 3])) ./ min(voxdim([2 3]));
    im2 = imresize(im2, [size(im2,1)*a size(im2,2)]);
end

size_max = max([max(size(im1)) max(size(im2)) max(size(im3))]);
three_plane_img = zeros([size_max 3*size_max]);
three_plane_img(:,1:size_max)              = image_center(im1, size_max);
three_plane_img(:,size_max*2+(1:size_max)) = image_center(im3, size_max);
three_plane_img(:,size_max+(1:size_max))   = image_center(im2, size_max);

MRS_struct.mask.img(ii,:,:)   = three_plane_img;
MRS_struct.mask.T1image(ii,:) = {nii_file};

warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:qhullmx:InternalWarning');

end



