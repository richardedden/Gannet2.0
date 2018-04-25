function MRS_struct = SiemensDICOMRead(MRS_struct,folder,waterfolder)
%% MRS_struct = SiemensDICOMRead(MRS_struct,folder,waterfolder)
%   This function is designed to load edited MR spectroscopy data in the 
%   Siemens flavour of DICOM data into a Gannet file structure. Files usually
%   have the extension '.IMA' or '.ima', and contain exactly 1 FID per
%   file, i.e. an acquisition of 320 averages will yield 320 IMA files.
%   
%   The user must specify the folder containing all of these averages. It
%   is assumed that they are ordered in the order of acquisition.
%   Water-suppressed and water-unsuppressed files need to be stored in
%   separate folders, which need to be defined accordingly:
%
%   Example:
%       folder = '/user/data/subject01/ima_gaba/';
%       waterfolder = '/user/data/subject01/ima_water/'; (optional)
%       MRS_struct = SiemensDICOMRead(MRS_struct,folder,waterfolder);
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2016-11-10)
%       goeltzs1@jhmi.edu
%   
%   Credits:    
% 
%   Version history:
%   0.9: First version (2016-11-10)
%   0.91: Added function for water data loading (2017-02-03)
%   0.92: Improved header parsing (2017-03-27). Thanks to Maria Yanez Lopez
%           and Ines Violante.
%   0.93: Added batch processing function (2017-11-16). Thanks to Dieter
%           Meyerhoff.
%   0.94: Added support for CMRR sequence (Eddie Auerbach, CMRR, University
%           of Minnesota) (2017-11-20). Thanks to Jim Lagopoulos.
%   0.95: Fills missing voxel geometry parameters in DICOM header with zero
%           values. Thanks to Alen Tersakyan.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% PREPARATION %%%
% Loop over number of datasets
ii = MRS_struct.ii;

% Locate folder and find all files in it. Will usually all be either upper or
% lower case, so concatenating the results of both should be fine and with
% no overlap.
% ima_file_list = [dir([folder,'/*.IMA']); dir([folder,'/*.ima'])]; % may
% cause problems on win/unix systems, take out for now % GO 11/16/2016
ima_file_list = dir([folder,'/*.IMA']); % GO 11/16/2016
fprintf('%d water-suppressed IMA files detected in %s.\n',length(ima_file_list),folder);

disp('Reading water-suppressed files...')

% Ordering of these files is not correct (i.e. 1,10,100,101...). Sort
% naturally.
ima_file_names = sort_nat({ima_file_list.name});
% Add folder to filenames (in case GannetLoad is run outside the folder)
% GO 11/20/2016
ima_file_names = strcat(folder, filesep, ima_file_names); % GO 11/20/2016
%%% /PREPARATION %%%

%%% HEADER INFO PARSING %%%
DicomHeader = read_dcm_header(ima_file_names{1});
MRS_struct.p.seq = DicomHeader.sequenceFileName;
MRS_struct.p.TR(ii) = DicomHeader.TR;
MRS_struct.p.TE(ii) = DicomHeader.TE;
MRS_struct.p.npoints(ii) = DicomHeader.vectorSize;
MRS_struct.p.Navg(ii) = 2*DicomHeader.nAverages;
MRS_struct.p.nrows(ii) = 2*DicomHeader.nAverages;
MRS_struct.p.sw(ii) = 1/DicomHeader.dwellTime * 1E9 * 0.5; % check with oversampling? hence factor 0.5, need to figure out why <=> probably dataset with 512 points, oversampled is 1024
MRS_struct.p.LarmorFreq(ii) = DicomHeader.tx_freq * 1E-6;
MRS_struct.p.voxdim(ii,1) = DicomHeader.VoI_PeFOV;
MRS_struct.p.voxdim(ii,2) = DicomHeader.VoI_RoFOV;
MRS_struct.p.voxdim(ii,3) = DicomHeader.VoIThickness;
MRS_struct.p.VoI_InPlaneRot(ii) = DicomHeader.VoI_InPlaneRot;
MRS_struct.p.voxoff(ii,1) = DicomHeader.PosSag;
MRS_struct.p.voxoff(ii,2) = DicomHeader.PosCor;
MRS_struct.p.voxoff(ii,3) = DicomHeader.PosTra;
MRS_struct.p.NormCor(ii) = DicomHeader.NormCor;
MRS_struct.p.NormSag(ii) = DicomHeader.NormSag;
MRS_struct.p.NormTra(ii) = DicomHeader.NormTra;
if isfield(MRS_struct.p, 'Siemens')
    MRS_struct.p.Siemens.deltaFreq = DicomHeader.deltaFreq;
    MRS_struct.p.Siemens.editRF.freq = DicomHeader.editRF.freq;
    MRS_struct.p.Siemens.editRF.bw = DicomHeader.editRF.bw;
    MRS_struct.p.Siemens.editRF.centerFreq = DicomHeader.editRF.centerFreq;
end
%%% /HEADER INFO PARSING %%%

%%% DATA LOADING %%%
% Preallocate array in which the FIDs are to be extracted.
MRS_struct.fids.data = zeros(MRS_struct.p.npoints(ii),length(ima_file_names));

% Collect all FIDs and sort them into MRS_struct
for kk = 1:length(ima_file_names)
    
    % Open IMA
    fd = dicom_open(ima_file_names{kk});
    
    % read the signal in as a complex FID
    MRS_struct.fids.data(:,kk) = dicom_get_spectrum_siemens(fd);

    fclose(fd);
end
disp('...complete')

% It appears that IMA stores the transients weirdly, 1-n/2 are all ONs, and
% n/2-n are all OFFS. Shuffle them below.
a = MRS_struct.fids.data(:,1:end/2);
b = MRS_struct.fids.data(:,1+end/2:end);
c = zeros(size(MRS_struct.fids.data));
c(:,1:2:end) = a;
c(:,2:2:end) = b;
MRS_struct.fids.data = c;
%%% /DATA LOADING %%%



%%% WATER DATA LOADING %%% % GO 02/05/2017
% If a water folder name is input to the function, repeat the same loading
% procedure for these files and hand the data over to the water data array
% of the MRS_struct.

% Set up the file name array.
if nargin == 3
    water_file_list = dir([waterfolder,'/*.IMA']);
    fprintf('%d water-unsuppressed IMA files detected in %s.\n',length(water_file_list),waterfolder);
    disp('Reading water-unsuppressed files...')
    water_file_names = sort_nat({water_file_list.name});
    water_file_names = strcat(waterfolder, filesep, water_file_names);
    
    % Load the actual water-unsuppressed data.
    MRS_struct.fids.waterdata = zeros(MRS_struct.p.npoints(ii),length(water_file_names));

    % Collect all FIDs and sort them into MRS_struct
    for kk = 1:length(water_file_names)
        
        % Open IMA
        fd = dicom_open(water_file_names{kk});
        
        % read the signal in as a complex FID
        MRS_struct.fids.data_water(:,kk) = dicom_get_spectrum_siemens(fd);
        
        fclose(fd);
    end
    disp('...complete')
end
%%% /WATER DATA LOADING %%%
