function MRS_struct = GannetMask_GE(fname, dcm_dir, MRS_struct, ii, vox, kk)

% MM (180329): Updated GannetMask_GE function. Rotated localizer files are
% no longer required. Voxel geometry is taken directly from P-file headers
% and the structural image DICOMs.
% Code heavily based on Ralph Noeske's (GE Berlin) SV_MRI voxel
% co-registration code.

% Parse P-file to extract voxel geometry
if MRS_struct.p.GE.rdbm_rev_num >= 11.0
    fid = fopen(fname, 'r', 'ieee-le');
else
    fid = fopen(fname, 'r', 'ieee-be');
end
fseek(fid, 1468, 'bof');
p_hdr_value = fread(fid, 12, 'integer*4');
fseek(fid, p_hdr_value(8), 'bof');
o_hdr_value = fread(fid, p_hdr_value(9)-p_hdr_value(8), 'real*4');
switch num2str(MRS_struct.p.GE.rdbm_rev_num)
    case '14.3'
        MRS_struct.p.voxdim(ii,:) = o_hdr_value(810:812)';
        MRS_struct.p.voxoff(ii,:) = o_hdr_value(813:815)';
        error('GannetMask_GE not yet compatible with rdbm_rev_num 14.3!');
    case '16'
        MRS_struct.p.voxdim(ii,:) = o_hdr_value(822:824)';
        MRS_struct.p.voxoff(ii,:) = o_hdr_value(825:827)';
        tlhc_LPS = o_hdr_value(905:907)';
        trhc_LPS = o_hdr_value(908:910)';
        brhc_LPS = o_hdr_value(911:913)';
    case '24'
        MRS_struct.p.voxdim(ii,:) = o_hdr_value(1228:1230)';
        MRS_struct.p.voxoff(ii,:) = o_hdr_value(1231:1233)';
        tlhc_LPS = o_hdr_value(1311:1313)';
        trhc_LPS = o_hdr_value(1314:1316)';
        brhc_LPS = o_hdr_value(1317:1319)';
    case '26.002'
        error('GannetMask_GE not yet compatible with rdbm_rev_num 26.002!');
end

MRS_struct.p.voxoff(ii,:) = MRS_struct.p.voxoff(ii,:) .* [-1 -1 1];
tlhc_LPS = tlhc_LPS .* [-1 -1 1];
trhc_LPS = trhc_LPS .* [-1 -1 1];
brhc_LPS = brhc_LPS .* [-1 -1 1];

e1_SVS_n = trhc_LPS - tlhc_LPS;
e1_SVS_n = e1_SVS_n ./ norm(e1_SVS_n);
e2_SVS_n = brhc_LPS - trhc_LPS;
e2_SVS_n = e2_SVS_n ./ norm(e2_SVS_n);
e3_SVS_n = -cross(e1_SVS_n, e2_SVS_n);
    
[~,orientation_SVS] = max(abs(e3_SVS_n));

if orientation_SVS == 3     % axial
    e1_SVS_n2 = e1_SVS_n;
    e2_SVS_n2 = e2_SVS_n;
    e3_SVS_n2 = e3_SVS_n;
elseif orientation_SVS == 2 % coronal
    e1_SVS_n2 = e1_SVS_n;
    e2_SVS_n2 = e3_SVS_n;
    e3_SVS_n2 = e2_SVS_n;
elseif orientation_SVS == 1 % sagittal
    e1_SVS_n2 = e3_SVS_n;
    e2_SVS_n2 = e1_SVS_n;
    e3_SVS_n2 = e2_SVS_n;
end

MRS_struct.p.voxang(ii,:) = get_euler(e1_SVS_n2, e2_SVS_n2, e3_SVS_n2);

e1_SVS = MRS_struct.p.voxdim(ii,1) * e1_SVS_n2;
e2_SVS = MRS_struct.p.voxdim(ii,2) * e2_SVS_n2;
e3_SVS = MRS_struct.p.voxdim(ii,3) * e3_SVS_n2;

% LPS gives center of voxel
LPS_SVS_edge = MRS_struct.p.voxoff(ii,:) - 0.5*e1_SVS - 0.5*e2_SVS - 0.5*e3_SVS;

% Read all DICOM files into one volume
currdir = pwd;
datadir = fileparts(MRS_struct.metabfile{ii});
if isempty(datadir)
    datadir = '.';
end
cd(datadir);
datadir = pwd;
cd(currdir);
cd(dcm_dir);
dcm_list = dir;
dcm_list = dcm_list(~ismember({dcm_list.name}, {'.','..','.DS_Store'}));
dcm_list = cellstr(char(dcm_list.name));
dcm_list = dcm_list(cellfun(@isempty, strfind(dcm_list, '.nii'))); %#ok<*STRCLFH>
dcm_list = dcm_list(cellfun(@isempty, strfind(dcm_list, '.mat')));
dcm_hdr = spm_dicom_headers(char(dcm_list));
nii_file_dir = spm_dicom_convert(dcm_hdr, 'all', 'flat', 'nii', datadir); % create NIFTI file of T1 image

% Rename NIFTI file into something more sensical (MM: 180810)
[a,~,c] = fileparts(nii_file_dir.files{1});
[~,b] = fileparts(MRS_struct.metabfile{ii});
nii_file = [a filesep b '_struc' c];
movefile(nii_file_dir.files{1}, nii_file);

slice_location = zeros(1,length(dcm_list));
for jj = 1:length(dcm_list)
    slice_location(jj) = dcm_hdr{jj}.SliceLocation;
end

% Order slices according to slice position
[~,order_index] = sort(slice_location);
tmp = dcm_hdr;
for jj = 1:length(dcm_list)
    dcm_hdr{jj} = tmp{order_index(jj)};    
end

cd(currdir);

V = spm_vol(nii_file);
T1 = spm_read_vols(V);

e1_MRI_n = dcm_hdr{1}.ImageOrientationPatient(1:3);
e2_MRI_n = dcm_hdr{1}.ImageOrientationPatient(4:6);
e3_MRI_n = cross(e1_MRI_n, e2_MRI_n); % e3 vector is perpendicular to the slice orientation

[~,orientation_MRI] = max(abs(e3_MRI_n));
if orientation_MRI == 2 % coronal
    e3_MRI_n = -e3_MRI_n;
end

MRI_voxel_size = [dcm_hdr{1}.PixelSpacing(1) ...
                  dcm_hdr{1}.PixelSpacing(2) ...
                  dcm_hdr{1}.SpacingBetweenSlices];

% LPS_edge gives location of the edge of the image volume
LPS_MRI_center = dcm_hdr{1}.ImagePositionPatient;
LPS_MRI_edge = LPS_MRI_center - 0.5 * MRI_voxel_size(1) * e1_MRI_n ...
                              - 0.5 * MRI_voxel_size(2) * e2_MRI_n ...
                              - 0.5 * MRI_voxel_size(3) * e3_MRI_n;

% Create voxel mask
E_MRI = [e1_MRI_n e2_MRI_n e3_MRI_n];
c_MRS = MRS_struct.p.voxoff(ii,:)';
c_MRI = E_MRI' * (c_MRS - LPS_MRI_edge);
d_MRI = c_MRI ./ MRI_voxel_size';
s_MRS = sqrt(sum(MRS_struct.p.voxdim(ii,:).^2))/2;
d_MRS = s_MRS ./ MRI_voxel_size';

[Xm,Ym,Zm] = ndgrid(1:V.dim(1), 1:V.dim(2), 1:V.dim(3));
X = LPS_MRI_center(1) + (Xm-1) * MRI_voxel_size(1) * e1_MRI_n(1) + (Ym-1) * MRI_voxel_size(2) * e2_MRI_n(1) + (Zm-1) * MRI_voxel_size(3) * e3_MRI_n(1);
Y = LPS_MRI_center(2) + (Xm-1) * MRI_voxel_size(1) * e1_MRI_n(2) + (Ym-1) * MRI_voxel_size(2) * e2_MRI_n(2) + (Zm-1) * MRI_voxel_size(3) * e3_MRI_n(2);
Z = LPS_MRI_center(3) + (Xm-1) * MRI_voxel_size(1) * e1_MRI_n(3) + (Ym-1) * MRI_voxel_size(2) * e2_MRI_n(3) + (Zm-1) * MRI_voxel_size(3) * e3_MRI_n(3);

P_1 = LPS_SVS_edge;
P_2 = LPS_SVS_edge + e1_SVS; % L
P_3 = LPS_SVS_edge + e2_SVS; % P
P_4 = LPS_SVS_edge + e3_SVS; % S
A = zeros(3,1);
mask = zeros(V.dim);

for e1 = max(floor(d_MRI(1) - d_MRS(1)), 0) : min(ceil(d_MRI(1) + d_MRS(1)), size(mask,1))         % L
    for e2 = max(floor(d_MRI(2) - d_MRS(2)), 0) : min(ceil(d_MRI(2) + d_MRS(2)), size(mask,2))     % P
        for e3 = max(floor(d_MRI(3) - d_MRS(3)), 0) : min(ceil(d_MRI(3) + d_MRS(3)), size(mask,3)) % S
            A(1) = X(e1,e2,e3);
            A(2) = Y(e1,e2,e3);
            A(3) = Z(e1,e2,e3);
            % Distance of A to planes in SI direction
            d_5 = e3_SVS_n2 * (A - P_1');
            d_6 = -e3_SVS_n2 * (A - P_4');
            if d_5 >= 0 && d_6 >= 0
                % Distance of A to planes in AP direction
                d_3 = e2_SVS_n2 * (A - P_1');
                d_4 = -e2_SVS_n2 * (A - P_3');
                if d_3 >= 0 && d_4 >= 0
                    % Distance of A to planes in RL direction
                    d_1 = e1_SVS_n2 * (A - P_1');
                    d_2 = -e1_SVS_n2 * (A - P_2');
                    if d_1 >= 0 && d_2 >= 0
                        mask(e1,e2,e3) = 1;
                    end
                end
            end
        end
    end
end

% LPS to RAS
mask = flip(mask,2);

% Output mask
[a,b] = fileparts(fname);
if isempty(a)
    a = '.';
end
V_mask.fname   = fullfile([a filesep b '_mask.nii']);
V_mask.descrip = 'MRS_voxel_mask';
V_mask.dim     = V.dim;
V_mask.dt      = V.dt;
V_mask.mat     = V.mat;
V_mask         = spm_write_vol(V_mask,mask);

MRS_struct.mask.(vox{kk}).outfile(ii,:) = cellstr(V_mask.fname);

% Transform structural image and co-registered voxel mask from voxel to
% world space for output (MM: 180221)
voxel_ctr = MRS_struct.p.voxoff(ii,:);
voxel_ctr(1:2) = -voxel_ctr(1:2);
[img_t,img_c,img_s] = voxel2world_space(V,voxel_ctr);
[mask_t,mask_c,mask_s] = voxel2world_space(V_mask,voxel_ctr);

img_t = flipud(img_t/max(T1(:)));
img_c = flipud(img_c/max(T1(:)));
img_s = flipud(img_s/max(T1(:)));

img_t = img_t + 0.175*flipud(mask_t);
img_c = img_c + 0.175*flipud(mask_c);
img_s = img_s + 0.175*flipud(mask_s);

size_max = max([max(size(img_t)) max(size(img_c)) max(size(img_s))]);
three_plane_img = zeros([size_max 3*size_max]);
three_plane_img(:,1:size_max)              = image_center(img_t, size_max);
three_plane_img(:,size_max+(1:size_max))   = image_center(img_s, size_max);
three_plane_img(:,size_max*2+(1:size_max)) = image_center(img_c, size_max);

MRS_struct.mask.(vox{kk}).img{ii} = three_plane_img;
MRS_struct.mask.(vox{kk}).T1image(ii,:) = {nii_file};

end


function euler_angles = get_euler(r1, r2, r3)

r1(3) = -r1(3);
r2(3) = -r2(3);
r3(3) = -r3(3);

if abs(r3(1)) ~= 1
    theta1 = -asin(r3(1));
    %theta2 = pi - theta1;
    psi1 = atan2(r3(2)/cos(theta1), r3(3)/cos(theta1));
    %psi2 = atan2(r3(2)/cos(theta2), r3(3)/cos(theta2));
    phi1 = atan2(r2(1)/cos(theta1), r1(1)/cos(theta1));
    %phi2 = atan2(r2(1)/cos(theta2), r1(1)/cos(theta2));
else
    phi1 = 0;
    if r3(1) == -1
        theta1 = pi/2;
        psi1 = phi1 + atan2(r1(2), r1(3));
    else
        theta1 = -pi/2;
        psi1 = -phi1 + atan2(-r1(2), -r1(3));
    end
end

euler_angles(1) = round(-phi1*180/pi);
euler_angles(2) = round(-psi1*180/pi);
euler_angles(3) = round(theta1*180/pi);

end



