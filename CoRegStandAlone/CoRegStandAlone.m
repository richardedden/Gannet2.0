function MRS_struct = CoRegStandAlone(metabfile,niifile)
% CoRegStandAlone(metabfile,niifile)
%   Reads the relevant header information from MRS files, performs 
%   SPM-based co-registration to a 3D anatomical image in nifti (*.nii)
%   format, and returns the tissue class fractions of gray matter, white
%   matter, and cerebrospinal fluid in the MRS voxel. 
%
%   Requires:
%       - SPM 12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
%
%   Input:
%       metabfile - cell containing the path to an MRS data file in one of 
%                   the following formats:
%                   Philips SDAT/SPAR, Siemens TWIX (*.dat), Siemens
%                   RDA (*.rda), GE (.7), DICOM (*.ima, *.dcm)
%       niifile   - cell containing the path to a *.nii file
%
%       It is possible to batch process data. In this case, the number of
%       elements of the metabfile cell needs to match the number of
%       elements of the niifile cell, with matching order.
%
%   Example:
%       MRS_struct = CoRegStandAlone({'MRS1.dat, 'MRS2.dat'}, {'nii.dat', 'nii.dat'})
%                                       
%       This example will co-register and segment two Siemens MRS voxels to
%       the same structural.
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-09-19)
%       goeltzs1@jhmi.edu
%   
%   History:
%       2018-09-19: First version of the code.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Pre-initialise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRS_struct.version.load = '180326'; % set to date when final updates have been made
MRS_struct.ii = 0;
MRS_struct.metabfile = metabfile;

% Flags
MRS_struct.p.mat = 1;       % Save results in *.mat output structure? (0 = NO, 1 = YES (default)).
MRS_struct.p.Vox = {'vox1'};  % Name of the voxel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2. Determine data parameters from header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discern input data format
MRS_struct = GannetDiscernDatatype(metabfile{1}, MRS_struct);

% Create output folder
if ~exist('CoRegStandAlone_output','dir')
    mkdir CoRegStandAlone_output;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   3. Load data from files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:length(metabfile) % Loop over all files in the batch (from metabfile)
    
    MRS_struct.ii = ii;
    
    switch MRS_struct.p.vendor
        
        case 'GE'
            MRS_struct = GERead(MRS_struct, metabfile{ii});
            
        case 'Siemens_twix'
            MRS_struct = SiemensTwixRead(MRS_struct, metabfile{ii});
            
        case 'Siemens_dicom' % GO 11/01/2016
            MRS_struct = SiemensDICOMRead(MRS_struct,metabfile{ii}); % GO 11/01/2016
            
        case 'dicom' % GO 11/30/2016
            MRS_struct = DICOMRead(MRS_struct,metabfile{ii}); % GO 11/01/2016
            
        case 'Siemens_rda'
            MRS_struct = SiemensRead(MRS_struct, metabfile{ii}, metabfile{ii});
            
        case 'Philips'
            MRS_struct = PhilipsRead(MRS_struct, metabfile{ii});
            
        case 'Philips_data'
            MRS_struct = PhilipsRead_data(MRS_struct, metabfile{ii});
        
        case 'Philips_raw'
            MRS_struct = PhilipsRawLoad(MRS_struct,metabfile{ii},3,0); % GO 11/02/2016 
    
    end % end of vendor switch loop for data load
        
end % end of load-and-processing loop over datasets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   4. Call coregister function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:length(metabfile) % Loop over all files in the batch (from metabfile)
    MRS_struct = CoReg(MRS_struct, niifile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   5. Call segment function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:length(metabfile) % Loop over all files in the batch (from metabfile)
    MRS_struct = Seg(MRS_struct);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   6. Clean up, save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Save MRS_struct as mat file
MRS_struct = rmfield(MRS_struct,'fids');
if ii == length(metabfile) && MRS_struct.p.mat
    % Set up filename
    mat_name = ['CoRegStandAlone_output/MRS_struct_CoRegStandAlone.mat'];
    save(mat_name,'MRS_struct');
end
           
end



