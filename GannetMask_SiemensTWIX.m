function MRS_struct = GannetMask_SiemensTWIX(filename, nii_file, MRS_struct, ii)

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:qhullmx:InternalWarning');

[path,name] = fileparts(filename);
fidoutmask = fullfile(path,[name '_mask.nii']);

% Extract voxel position and rotation parameters from MRS_struct
NormSag = MRS_struct.p.NormSag(ii);
NormCor = MRS_struct.p.NormCor(ii);
NormTra = MRS_struct.p.NormTra(ii);
VoxOffs = [MRS_struct.p.voxoff(ii,1) MRS_struct.p.voxoff(ii,2) MRS_struct.p.voxoff(ii,3)];
VoxDims = [MRS_struct.p.voxdim(ii,1) MRS_struct.p.voxdim(ii,2) MRS_struct.p.voxdim(ii,3)];
VoI_InPlaneRot = MRS_struct.p.VoI_InPlaneRot(ii);

% Define normal vector, rotation angle, and orientation of the imaging slab
Norm = [-NormSag -NormCor NormTra];
%Norm = Norm *-1;
ROT = VoI_InPlaneRot;
sl_orient = 's'; % 't' = transversal, 's' = sagittal', 'c' = coronal;

% phase reference vector - 
% adapted from Andre van der Kouwe's "autoaligncorrect.cpp"
Phase	= zeros(3, 1);
switch sl_orient
    case 't'
        Phase(1)	= 0;
        Phase(2)	=  Norm(3)*sqrt(1/(Norm(2)*Norm(2)+Norm(3)*Norm(3)));
        Phase(3)	= -Norm(2)*sqrt(1/(Norm(2)*Norm(2)+Norm(3)*Norm(3)));
    case 'c'
        Phase(1)	=  Norm(2)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(2)	= -Norm(1)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(3)	= 0;
    case 's'
        Phase(1)	= -Norm(2)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(2)	=  Norm(1)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(3)	= 0;
    otherwise
        fprintf(1, 'Unknown orientation parameter passed. Returning with dummy M_R');
	return;
end

% The readout reference vector is the cross product of Norm and Phase
Readout = cross(Norm, Phase);
M_R = zeros(4, 4);
M_R(1:3, 1)	= Phase;
M_R(1:3, 2)	= Readout;
M_R(1:3, 3) = Norm;

% Define matrix for rotation around in-plane rotation angle
M3_Mu	= [	 cos(ROT)	sin(ROT)	0
            -sin(ROT)	cos(ROT)	0
            0           0           1];
        
M3_R	= M_R(1:3,1:3)	* M3_Mu;
M_R(1:3,1:3)	= M3_R;

% The MGH vox2ras matrix inverts the Readout column
M_R		= M_R *   [ 1  0  0  0
                    0 -1  0  0
                    0  0  1  0
                    0  0  0  1];

% Final rotation matrix
rotmat = M_R(1:3,1:3);

V = spm_vol(nii_file);
[T1,XYZ] = spm_read_vols(V);

%Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
%tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
[~,voxdim] = spm_get_bbox(V,'fv'); % MM (180220)
voxdim = abs(voxdim)';
halfpixshift = -voxdim(1:3)/2;
halfpixshift(3) = -halfpixshift(3);
XYZ = XYZ + repmat(halfpixshift, [1 size(XYZ,2)]);

% We need to flip ap and lr axes to match NIFTI convention
VoxOffs(1) = -VoxOffs(1);
VoxOffs(2) = -VoxOffs(2);

% Define voxel coordinates before rotation and transition
vox_ctr = ...
    [VoxDims(1)/2 -VoxDims(2)/2  VoxDims(3)/2;
    -VoxDims(1)/2 -VoxDims(2)/2  VoxDims(3)/2;
    -VoxDims(1)/2  VoxDims(2)/2  VoxDims(3)/2;
     VoxDims(1)/2  VoxDims(2)/2  VoxDims(3)/2;
    -VoxDims(1)/2  VoxDims(2)/2 -VoxDims(3)/2;
     VoxDims(1)/2  VoxDims(2)/2 -VoxDims(3)/2;
     VoxDims(1)/2 -VoxDims(2)/2 -VoxDims(3)/2;
    -VoxDims(1)/2 -VoxDims(2)/2 -VoxDims(3)/2];

% Apply rotation as prescribed
vox_rot = rotmat*vox_ctr.';

% Shift rotated voxel by the center offset to its final position
vox_ctr_coor = [VoxOffs(1) VoxOffs(2) VoxOffs(3)];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
vox_corner = vox_rot + vox_ctr_coor;

% Create a mask with all voxels that are inside the voxel
mask = zeros(1,size(XYZ,2));
sphere_radius = sqrt((VoxDims(1)/2)^2+(VoxDims(2)/2)^2+(VoxDims(3)/2)^2);
distance2voxctr = sqrt(sum((XYZ-repmat([VoxOffs(1) VoxOffs(2) VoxOffs(3)].',[1 size(XYZ, 2)])).^2,1));
sphere_mask(distance2voxctr <= sphere_radius) = 1;
mask(sphere_mask == 1) = 1;
XYZ_sphere = XYZ(:,sphere_mask == 1);
tri = delaunayn([vox_corner.'; [VoxOffs(1) VoxOffs(2) VoxOffs(3)]]);
tn = tsearchn([vox_corner.'; [VoxOffs(1) VoxOffs(2) VoxOffs(3)]], tri, XYZ_sphere.');
isinside = ~isnan(tn);
mask(sphere_mask==1) = isinside;

% Take over the voxel dimensions from the structural
mask = reshape(mask, V.dim);

V_mask.fname   = fidoutmask ;
V_mask.descrip = 'MRS_voxel_mask';
V_mask.dim     = V.dim;
V_mask.dt      = V.dt;
V_mask.mat     = V.mat;

spm_write_vol(V_mask,mask);

% Build output page
fidoutmask = cellstr(fidoutmask);
MRS_struct.mask.outfile(ii,:) = fidoutmask;
% Not clear how to formulate the rotations for triple rotations (revisit)
MRS_struct.p.voxang(ii,:) = [NaN NaN NaN];  

% Transform structural image and co-registered voxel mask from voxel to
% world space for output (MM: 180221)
[img_t,img_c,img_s] = voxel2world_space(V,VoxOffs);
V_mask = spm_vol(fidoutmask{1});
[mask_t,mask_c,mask_s] = voxel2world_space(V_mask,VoxOffs);

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

MRS_struct.mask.img{ii} = three_plane_img;
MRS_struct.mask.T1image(ii,:) = {nii_file};

warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:qhullmx:InternalWarning');

end



