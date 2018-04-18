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
% Simply open the first IMA file, the information should all be the same.
fid = fopen(ima_file_names{1});

% Start looking for a convenient parameter block in the first IMA file. The
% line before will start with ### ASCCONV BEGIN and end with ### ASCCONV
% END. This is defined here.
head_start_text = '### ASCCONV BEGIN';
head_end_text   = '### ASCCONV END';

tline = fgets(fid);

while (isempty(strfind(tline , head_end_text))) %#ok<*STREMP>
    
    tline = fgets(fid);
    
    if ( isempty(strfind (tline , head_start_text)) + isempty(strfind (tline , head_end_text )) == 2)
                
        % Find lines with 'equal' signs, all the information is in there.    
        findequal = strfind(tline,'=');
        variable = strtrim(tline(1:findequal-1)) ;
        value    = strtrim(tline(findequal+1 : length(tline))) ;
        
        switch variable
            case {'tSequenceFileName'}
                MRS_struct.p.seq = value;
            case {'alTR[0]'}
                MRS_struct.p.TR(ii) = str2double(value) / 1000;
            case {'alTE[0]'}
                MRS_struct.p.TE(ii) = str2double(value) / 1000;
            case {'sSpecPara.lVectorSize'}
                MRS_struct.p.npoints(ii) = str2double(value);
            case {'lAverages'}
                % Minnesota sequence (CMRR, Eddy Auerbach) stores numbers of averages in a
                % different field. GO 112017.
                if strcmp(MRS_struct.p.seq,'""%CustomerSeq%\eja_svs_mpress""')
                    % Don't bother
                else
                    MRS_struct.p.Navg(ii) = 2*str2double(value);
                    MRS_struct.p.nrows(ii) = 2*str2double(value);
                end
            case {'sWipMemBlock.alFree[2]'}
                % Minnesota sequence (CMRR, Eddy Auerbach) stores numbers of averages in a
                % different field. GO 112017.
                if strcmp(MRS_struct.p.seq,'""%CustomerSeq%\eja_svs_mpress""')
                    MRS_struct.p.Navg(ii) = 2*str2double(value);
                    MRS_struct.p.nrows(ii) = 2*str2double(value);
                else
                    % Don't bother
                end
            case {'sRXSPEC.alDwellTime[0]'}
                MRS_struct.p.sw(ii) = (1/str2double(value)) * 1E9 * 0.5; % check with oversampling? hence factor 0.5, need to figure out why <=> probably dataset with 512 points, oversampled is 1024
            case {'lFrequency','sTXSPEC.asNucleusInfo[0].lFrequency'}
                MRS_struct.p.LarmorFreq(ii) = str2double(value) * 1E-6;
            case {'sSpecPara.sVoI.dPhaseFOV'}
                MRS_struct.p.voxdim(ii,1) = str2double(value);
            case {'sSpecPara.sVoI.dReadoutFOV'}
                MRS_struct.p.voxdim(ii,2) = str2double(value);
            case {'sSpecPara.sVoI.dThickness'}
                MRS_struct.p.voxdim(ii,3) = str2double(value);
            case {'sSpecPara.sVoI.dInPlaneRot'}
                MRS_struct.p.VoI_InPlaneRot(ii) = str2double(value);
            case {'sSpecPara.sVoI.sPosition.dSag'}
                MRS_struct.p.voxoff(ii,1) = str2double(value);
            case {'sSpecPara.sVoI.sPosition.dCor'}
                MRS_struct.p.voxoff(ii,2) = str2double(value);
            case {'sSpecPara.sVoI.sPosition.dTra'}
                MRS_struct.p.voxoff(ii,3) = str2double(value);
            case {'sSpecPara.sVoI.sNormal.dCor'}
                MRS_struct.p.NormCor(ii) = str2double(value);
            case {'sSpecPara.sVoI.sNormal.dSag'}
                MRS_struct.p.NormSag(ii) = str2double(value);
            case {'sSpecPara.sVoI.sNormal.dTra'}
                MRS_struct.p.NormTra(ii) = str2double(value);
            case {'sSpecPara.dDeltaFrequency'}
                MRS_struct.p.Siemens.deltaFreq = str2double(value);
            case {'sWipMemBlock.adFree[7]'}    
                MRS_struct.p.Siemens.editRF.freq = str2double(value);
            case {'sWipMemBlock.adFree[8]'}
                MRS_struct.p.Siemens.editRF.bw = str2double(value);
            case {'sWipMemBlock.adFree[9]'}
                MRS_struct.p.Siemens.editRF.centerFreq = str2double(value);
            otherwise
                % don't care, do nothing
        end
    else
        % don't care about the rest, do nothing
    end
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
