function MRS_struct = GannetMask_SiemensRDA(filename, nii_file, MRS_struct, ii, vox, kk)

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:qhullmx:InternalWarning');

% Still kind of beta with updates from GO 2017
% this relies on SPM and a nifti file
% being testing on data from Univ of Florida - will need to extend
% also may need to extend to change dicoms into nifti

if nargin == 2
    MRS_struct.ii = 1;
    ii = 1;
end

% Parse RDA filename and establish nifti voxelmask filename
[path,name,~] = fileparts(filename);
fidoutmask = fullfile(path,[name '_mask.nii']);
fid = fopen(filename);
disp(filename);

% Go through RDA header line by line and extract header info
head_start_text = '>>> Begin of header <<<';
head_end_text   = '>>> End of header <<<';
tline = fgets(fid);

while (isempty(strfind(tline, head_end_text))) %#ok<*STREMP>
    
    tline = fgets(fid);
    
    if ( isempty(strfind(tline, head_start_text)) + isempty(strfind(tline, head_end_text)) == 2)
                
        % Store this data in the appropriate format
        occurence_of_colon = strfind(tline,':');
        variable = tline(1:occurence_of_colon-1);
        value    = tline(occurence_of_colon+1 : length(tline)); 
        
        switch variable
        case { 'VOINormalSag' , 'VOINormalCor' , 'VOINormalTra' , 'VOIPositionSag', 'VOIPositionCor', 'VOIPositionTra', 'VOIThickness','VOIReadoutFOV','VOIPhaseFOV'  }
            eval(['rda.' , variable , ' = str2num(value); ']);
        case { 'RowVector[0]' }
            rda.row(1)=str2double(value);
        case { 'RowVector[1]' }
            rda.row(2)=str2double(value);
        case { 'RowVector[2]' }
            rda.row(3)=str2double(value);
        case { 'ColumnVector[0]' }
            rda.column(1)=str2double(value);
        case { 'ColumnVector[1]' }
            rda.column(2)=str2double(value);
        case { 'ColumnVector[2]' }
            rda.column(3)=str2double(value);
        case { 'PositionVector[0]' }
            rda.position(1)=str2double(value);
        case { 'PositionVector[1]' }
            rda.position(2)=str2double(value);
        case { 'PositionVector[2]' }
            rda.position(3)=str2double(value);                  
        case {'VOIPositionSag' }
            rda.pdSag = str2double(value);
        case {'VOIPositionCor' }
            rda.pdCor = str2double(value);
        case {'VOIPositionTra' }
            rda.pdTra = str2double(value);
        end
        
    else
        % Don't bother storing this bit of the output
    end
    
    
end

fclose(fid);

% Create voxel coordinates
MRS_struct.p.voxoff(ii) =[ rda.VOIPositionSag rda.VOIPositionCor rda.VOIPositionTra];
MRS_struct.p.voxdim(ii) = [rda.VOIThickness rda.VOIReadoutFOV rda.VOIPhaseFOV ]; 
MRS_Rot(:,1)=rda.row.'.* [-1 -1 1]' ;
MRS_Rot(:,2)=rda.column.' ;
MRS_Rot(:,3)=cross(MRS_Rot(:,1),MRS_Rot(:,2));
MRS_Rot(1,:)=-MRS_Rot(1,:);
rotmat=-MRS_Rot;

% Change to current working directory
currdir=pwd;
cd(currdir);

% Open nifti file
V = spm_vol(nii_file);
[T1,XYZ] = spm_read_vols(V);

%Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
%tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
[~,voxdim] = spm_get_bbox(V,'fv'); % MM (180220)
voxdim = abs(voxdim)';
halfpixshift = -voxdim(1:3)/2;
halfpixshift(3) = -halfpixshift(3);
XYZ = XYZ+repmat(halfpixshift,[1 size(XYZ,2)]);

% Parse voxel dimensions
ap_size = MRS_struct.p.voxsize(2);
lr_size = MRS_struct.p.voxsize(1);
cc_size = MRS_struct.p.voxsize(3);
ap_off = MRS_struct.p.voxoff(2);
lr_off = MRS_struct.p.voxoff(1);
cc_off = MRS_struct.p.voxoff(3);
%ap_ang = MRS_struct.p.voxang(2);
%lr_ang = MRS_struct.p.voxang(1);
%cc_ang = MRS_struct.p.voxang(3);

%We need to flip ap and lr axes to match NIFTI convention
ap_off = -ap_off;
lr_off = -lr_off;

% Rotation may be required backwards, doublecheck
%ap_ang = -ap_ang;
%lr_ang = -lr_ang;

% define the voxel - use x y z  
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
   
vox_rot=rotmat*vox_ctr.';

% calculate corner coordinates relative to xyz origin
vox_ctr_coor = [lr_off ap_off cc_off];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
vox_corner = vox_rot+vox_ctr_coor;

% create mask
mask = zeros(1,size(XYZ,2));
sphere_radius = sqrt((lr_size/2)^2+(ap_size/2)^2+(cc_size/2)^2);
distance2voxctr=sqrt(sum((XYZ-repmat([lr_off ap_off cc_off].',[1 size(XYZ, 2)])).^2,1));
sphere_mask(distance2voxctr<=sphere_radius)=1;
mask(sphere_mask==1) = 1;
XYZ_sphere = XYZ(:,sphere_mask == 1);
tri = delaunayn([vox_corner.'; [lr_off ap_off cc_off]]);
tn = tsearchn([vox_corner.'; [lr_off ap_off cc_off]], tri, XYZ_sphere.');
isinside = ~isnan(tn);
mask(sphere_mask==1) = isinside;
mask = reshape(mask, V.dim);

% Fill mask header
V_mask.fname = fidoutmask;
V_mask.descrip = 'MRS_voxel_mask';
V_mask.dim = V.dim;
V_mask.dt = V.dt;
V_mask.mat = V.mat;

% Write to file
V_mask = spm_write_vol(V_mask,mask);

% construct output 
voxel_ctr = [-lr_off -ap_off cc_off];

% Populate MRS_struct with mask information
fidoutmask = cellstr(fidoutmask);
MRS_struct.mask.(vox{kk}).outfile(MRS_struct.ii,:) = fidoutmask;
MRS_struct.p.voxang(ii,:) = [NaN NaN NaN];  % put as NaN for now - for output page

voxel_ctr(1:2) = -voxel_ctr(1:2);

% Transform structural image and co-registered voxel mask from voxel to
% world space for output (MM: 180221)
[img_t,img_c,img_s] = voxel2world_space(V,voxel_ctr);
[mask_t,mask_c,mask_s] = voxel2world_space(V_mask,voxel_ctr);

img_t = flipud(img_t/max(T1(:)));
img_c = flipud(img_c/max(T1(:)));
img_s = flipud(img_s/max(T1(:)));

img_t = img_t + 0.08*flipud(mask_t);
img_c = img_c + 0.08*flipud(mask_c);
img_s = img_s + 0.08*flipud(mask_s);

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



