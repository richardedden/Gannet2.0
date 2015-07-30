function [MRS_struct ] = GannetMask_SiemensTWIX(filename, nii_file, MRS_struct, ii)

% this relies on SPM and a nifti t
% being testing on data from Univ of Florida - will need to extend
% also may need to extend to change dicoms into nifti


if(nargin == 2)
    MRS_struct.ii=1;
    ii = 1;
end
% some code to make ok wth 2 or 3  inputs
% fix nifti inputs and writing output mask to directory of where the nifti
% is. 

[path,name,ext] = fileparts(filename);
[pathnii,namenii,extnii] = fileparts(nii_file);

fidoutmask = fullfile(path,[name '_mask.nii'])

fid=fopen(filename);
line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_Normal_Sag">  { <Precision> ');
equals_index=strfind(line,'16 ');
            while isempty(index) || isempty(equals_index)
                line=fgets(fid);
                index=strfind(line,'<ParamDouble."VoI_Normal_Sag">  { <Precision> ');
                equals_index=strfind(line,'16 ');
            end
NormSag = line(equals_index+4:equals_index+16);
NormSag = str2double(NormSag);


line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_Normal_Cor">  { <Precision> ');
equals_index=strfind(line,'16 ');
            while isempty(index) || isempty(equals_index)
                line=fgets(fid);
                index=strfind(line,'<ParamDouble."VoI_Normal_Cor">  { <Precision> ');
                equals_index=strfind(line,'16 ');
            end
NormCor = line(equals_index+4:equals_index+16);
NormCor = str2double(NormCor);

line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_Normal_Tra">  { <Precision> ');
equals_index=strfind(line,'16 ');
            while isempty(index) || isempty(equals_index)
                line=fgets(fid);
                index=strfind(line,'<ParamDouble."VoI_Normal_Tra">  { <Precision> ');
                equals_index=strfind(line,'16 ');
            end
NormTra = line(equals_index+4:equals_index+16);
NormTra = str2double(NormTra);


line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_Position_Sag">  { <Precision> ');
equals_index=strfind(line,'16 ');
            while isempty(index) || isempty(equals_index)
                line=fgets(fid);
                index=strfind(line,'<ParamDouble."VoI_Position_Sag">  { <Precision> ');
                equals_index=strfind(line,'16 ');
            end
PosSag = line(equals_index+4:equals_index+16);
PosSag = str2double(PosSag);


line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_Position_Cor">  { <Precision> ');
equals_index=strfind(line,'16 ');
            while isempty(index) || isempty(equals_index)
                line=fgets(fid);
                index=strfind(line,'<ParamDouble."VoI_Position_Cor">  { <Precision> ');
                equals_index=strfind(line,'16 ');
            end
PosCor = line(equals_index+4:equals_index+16);
PosCor = str2double(PosCor);

line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_Position_Tra">  { <Precision> ');
equals_index=strfind(line,'16 ');
            while isempty(index) || isempty(equals_index)
                line=fgets(fid);
                index=strfind(line,'<ParamDouble."VoI_Position_Tra">  { <Precision> ');
                equals_index=strfind(line,'16 ');
            end
PosTra = line(equals_index+4:equals_index+16);
PosTra = str2double(PosTra);

line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_SliceThickness">  { <Precision> ');
equals_index=strfind(line,'16 ');
            while isempty(index) || isempty(equals_index)
                line=fgets(fid);
                index=strfind(line,'<ParamDouble."VoI_SliceThickness">  { <Precision> ');
                equals_index=strfind(line,'16 ');
            end
VOIThickness = line(equals_index+4:equals_index+16);
VOIThickness = str2double(VOIThickness);

fclose(fid);
fid=fopen(filename);
line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_InPlaneRotAngle">  { <Precision> ');
equals_index=strfind(line,'16 ');
            while isempty(index) || isempty(equals_index)
                line=fgets(fid);
                index=strfind(line,'<ParamDouble."VoI_InPlaneRotAngle">  { <Precision> ');
                equals_index=strfind(line,'16 ');
            end
VoI_InPlaneRot = line(equals_index+4:equals_index+16);
VoI_InPlaneRot = str2double(VoI_InPlaneRot);


line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_RoFOV">  { <Precision> ');
equals_index=strfind(line,'16 ');
            while isempty(index) || isempty(equals_index)
                line=fgets(fid);
                index=strfind(line,'<ParamDouble."VoI_RoFOV">  { <Precision> ');
                equals_index=strfind(line,'16 ');
            end
VoI_RoFOV = line(equals_index+4:equals_index+16);
VoI_RoFOV = str2double(VoI_RoFOV);

line = fgets(fid);
index = strfind(line, '<ParamDouble."VoI_PeFOV">  { <Precision> ');
equals_index=strfind(line,'16 ');
            while isempty(index) || isempty(equals_index)
                line=fgets(fid);
                index=strfind(line,'<ParamDouble."VoI_PeFOV">  { <Precision> ');
                equals_index=strfind(line,'16 ');
            end
VoI_PeFOV = line(equals_index+4:equals_index+16);
VoI_PeFOV = str2double(VoI_PeFOV);


fclose(fid);

%%

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
Raxisangle(4,:) = [ 0 0 0 1];

colStart(1,1)=0;
colStart(3,1)=1/sqrt(((-1.0 * ZED(3) / ZED(2))^2)+1);
colStart(2,1)=-1.0 * colStart(3) * ZED(3) / ZED(2);

ColStart = colStart;
ColStart(1:3,4) = 0;
ColStart(4,:) = [0 0 0 1];

ColMat = Raxisangle*ColStart;
Col(1) = ColMat(1,1);
Col(2) = ColMat(2,1);
Col(3) = ColMat(3,1);

Col

Row(1) = Col(2) * ZED(3) - Col(3) * ZED(2);
Row(2) = Col(3) * ZED(1) - Col(1) * ZED(3);
Row(3) = Col(1) * ZED(2) - Col(2) * ZED(1);

Row



%MRS_struct.p.voxoff=[ rda.position(1) rda.position(2) rda.position(3)];
MRS_struct.p.voxoff=[PosSag PosCor PosTra];
MRS_struct.p.voxsize = [VoI_RoFOV VoI_PeFOV VOIThickness]; 
MRS_Rot(:,1) = Row';
MRS_Rot(:,2) = Col';

MRS_Rot(:,1)= MRS_Rot(:,1) .* [ -1 -1 1]'
MRS_Rot(:,2)= MRS_Rot(:,2) .* [ -1 -1 1]'

MRS_Rot(:,3)=cross(MRS_Rot(:,2),MRS_Rot(:,1))

rotmat=MRS_Rot % used to be rotmat = -MRS_rot but doesn't seem to do anything 
    %- bc of corners defn in mat - the values just get reorded in the vox_corner mat
%rotmat(:,3) = [1 1 1];
%error('gothere');

% currdir=pwd;
% cd(currdir);

%voxang is initialized to zero so the code runs, but ideally, it needs to
%parse the rotation info from nii_file2.



V=spm_vol(nii_file);
%V.dim
[T1,XYZ]=spm_read_vols(V);
H=spm_read_hdr(nii_file);

%Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
%tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
halfpixshift = -H.dime.pixdim(1:3).'/2;
halfpixshift(3) = -halfpixshift(3);
XYZ=XYZ+repmat(halfpixshift,[1 size(XYZ,2)]);

ap_size = MRS_struct.p.voxsize(2);
lr_size = MRS_struct.p.voxsize(1);
cc_size = MRS_struct.p.voxsize(3);
ap_off = MRS_struct.p.voxoff(2);
lr_off = MRS_struct.p.voxoff(1);
cc_off = MRS_struct.p.voxoff(3);
%ap_ang = MRS_struct.p.voxang(2);
%lr_ang = MRS_struct.p.voxang(1);
%cc_ang = MRS_struct.p.voxang(3);
% 
% 
%We need to flip ap and lr axes to match NIFTI convention
ap_off = -ap_off;
lr_off = -lr_off;

%LINE 118

%ap_ang = -ap_ang;
%lr_ang = -lr_ang;



% define the voxel - use x y z  
% currently have spar convention - will need to
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

% vox_ctr = ...
%       [lr_size 0 cc_size ;
%       0 0 cc_size ;
%        0 ap_size cc_size ;
%        lr_size ap_size cc_size ;
%        0 ap_size 0 ;
%        lr_size ap_size 0 ;
%        lr_size 0 0 ;
%        0 0 0];

   
vox_rot=rotmat*vox_ctr.'

% calculate corner coordinates relative to xyz origin
vox_ctr_coor = [lr_off ap_off cc_off];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
vox_corner = vox_rot+vox_ctr_coor



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

%FOR NOW NEED TO FIX
%MRS_struct.mask.outfile(MRS_struct.ii,:)=fidoutmask;

voxel_ctr(1:2)=-voxel_ctr(1:2);
voxel_search=(XYZ(:,:)-repmat(voxel_ctr.',[1 size(XYZ,2)])).^2;
voxel_search=sqrt(sum(voxel_search,1));
[min2,index1]=min(voxel_search);
[slice(1) slice(2) slice(3)]=ind2sub( V.dim,index1);

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

MRS_struct.mask.img(MRS_struct.ii,:,:)=three_plane_img;


figure(198)
imagesc(three_plane_img);
colormap('gray');
caxis([0 1])
axis equal;
axis tight;
axis off;

%end

