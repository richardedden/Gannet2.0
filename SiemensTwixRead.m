function MRS_struct = SiemensTwixRead(MRS_struct,fname,fname_water)
%% function [ MRS_struct ] = SiemensTwixRead(MRS_struct,fname,fname_water)
%   Reads Siemens TWIX files (*.dat).
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2017-03-22)
%       goeltzs1@jhmi.edu
%   
%   Credits:    
%
%   History:
%       2017-03-22: First version.
%       2017-04-21: Move loading module to separate function, add 
%       support for loading PRESS water reference data.
%       2017-07-13: - Metabolite spectra phased according to unsuppressed
%                     MEGA-PRESS water reference acquisition
%                   - Make parsing of editing pulse frequencies available
%                     only when the fields are actually present (may depend 
%                     on vendor and sequence version). 
%                   - Minor improvements.

ii = MRS_struct.ii;

% Get the raw data and header info from the MEGA-PRESS files.
[MetabData, MetabHeader] = GetTwixData(fname);
% Populate MRS_struct with relevant info.
MRS_struct.p.pointsBeforeEcho       = MetabHeader.pointsBeforeEcho;
MRS_struct.p.sw                     = 1/MetabHeader.dwellTime;
MRS_struct.p.LarmorFreq             = MetabHeader.Bo*42.577;
MRS_struct.p.TR(ii)                 = MetabHeader.TR;
MRS_struct.p.TE(ii)                 = MetabHeader.TE;
MRS_struct.p.npoints                = size(MetabData,2);
MRS_struct.p.nrows                  = size(MetabData,3);
MRS_struct.p.Navg(ii)               = size(MetabData,3);
if exist('MetabHeader.editRF','var')
    MRS_struct.p.Siemens.editRF.freq(ii,:)      = TwixHeader.editRF.freq;
    MRS_struct.p.Siemens.editRF.centerFreq(ii)  = TwixHeader.editRF.centerFreq;
    MRS_struct.p.Siemens.editRF.bw(ii)          = TwixHeader.editRF.bw;
end
if exist('MetabHeader.deltaFreq','var')
    MRS_struct.p.Siemens.deltaFreq.metab(ii)    = TwixHeader.deltaFreq;
    MRS_struct.p.Siemens = reorderstructure(MRS_struct.p.Siemens, 'editRF', 'deltaFreq');
end

% If additional data points have been acquired before the echo starts,
% remove these here.
MetabData = MetabData(:,(MRS_struct.p.pointsBeforeEcho+1):end,:);
MRS_struct.p.npoints = MRS_struct.p.npoints - MRS_struct.p.pointsBeforeEcho; % MM (160914)

% If water reference is provided, load this one as well, and populate
% MRS_struct with water reference specific information.
if nargin == 3
    [WaterData, WaterHeader] = GetTwixData(fname_water);
    MRS_struct.p.pointsBeforeEcho_water  = WaterHeader.pointsBeforeEcho;
    MRS_struct.p.sw_water                = 1/WaterHeader.dwellTime;
    MRS_struct.p.TR_water(ii)            = WaterHeader.TR;
    MRS_struct.p.TE_water(ii)            = WaterHeader.TE;
    MRS_struct.p.npoints_water           = size(WaterData,2);
    MRS_struct.p.nrows_water             = size(WaterData,3);
    MRS_struct.p.Nwateravg(ii)           = size(WaterData,3);
    MRS_struct.p.seqtype_water           = WaterHeader.seqtype;
    if exist('WaterHeader.deltaFreq','var')
        MRS_struct.p.Siemens.deltaFreq.water(ii)    = TwixHeader.deltaFreq;
        MRS_struct.p.Siemens = reorderstructure(MRS_struct.p.Siemens, 'editRF', 'deltaFreq');
    end
    
    % If additional data points have been acquired before the echo starts,
    % remove these here.
    WaterData = WaterData(:,(MRS_struct.p.pointsBeforeEcho_water+1):end,:);
    MRS_struct.p.npoints_water = MRS_struct.p.npoints_water - MRS_struct.p.pointsBeforeEcho_water; % MM (160914)
    
    %Combine data based upon first point of FIDs (mean over all averages)
    % MM (170123)
    firstpoint_water = conj(WaterData(:,1,:));
    channels_scale = squeeze(sqrt(sum(firstpoint_water .* conj(firstpoint_water),1)));
    channels_scale = repmat(channels_scale, [1 size(WaterData,1) MRS_struct.p.npoints_water]);
    channels_scale = permute(channels_scale, [2 3 1]);
    firstpoint_water = repmat(firstpoint_water, [1 MRS_struct.p.npoints_water 1])./channels_scale;
    WaterData = WaterData .* firstpoint_water;
    WaterData = conj(squeeze(sum(WaterData,1)));
    WaterData = squeeze(mean(WaterData,2));
    MRS_struct.fids.data_water = double(WaterData);
end

% Phasing of metabolite data.
% If the water reference has been acquired with MEGA-PRESS, use the phase
% information from it to phase the metabolite data.
if isfield(MRS_struct.p,'seqtype_water') && strcmp(MRS_struct.p.seqtype_water,'MEGAPRESS')
    disp('MEGA-PRESS water reference found!');
    disp('Phasing metabolite data with water reference phase...');
    % Use first point of water data to phase water-suppressed data
    firstpoint = mean(firstpoint_water,3);
    firstpoint = repmat(firstpoint, [1 1 size(MetabData,3)]);
    MetabData = MetabData .* firstpoint;
    MetabData = conj(squeeze(sum(MetabData,1)));
    MRS_struct.fids.data = MetabData;
else
    % If no water data (or PRESS water reference) provided, combine data 
    % based upon first point of metabolite data (average all transients)
    if isfield(MRS_struct.p,'seqtype_water') && strcmp(MRS_struct.p.seqtype_water,'PRESS')
        disp('PRESS water reference found!');
    else
        disp('No water reference found!');
    end
    disp('Phasing metabolite data...');
    %Combine data based upon first point of FIDs (mean over all averages)
    firstpoint=mean(conj(MetabData(:,1,:)),3);
    channels_scale=squeeze(sqrt(sum(firstpoint.*conj(firstpoint))));
    firstpoint=repmat(firstpoint, [1 MRS_struct.p.npoints MRS_struct.p.nrows])/channels_scale;
    % Multiply the Multichannel data by the firstpointvector
    % zeroth order phasing of spectra
    MetabData = MetabData .* firstpoint;
    % sum over Rx channels
    MetabData = conj(squeeze(sum(MetabData,1)));
    MRS_struct.fids.data = double(MetabData);
end
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SEPARATE FUNCTIONS START BELOW %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TwixData, TwixHeader] = GetTwixData(fname)

% Pull TWIX data in with the mapVBVD tool
twix_obj=mapVBVD(fname);
            
% Is the data single-RAID or multi-RAID?
% struct - single-RAID
% cell - multi-RAID, with info in the last cell element
if isstruct(twix_obj)
    disp('loading single-RAID file...')
elseif iscell(twix_obj)
    disp('loading multi-RAID file...')
    twix_obj = twix_obj{end};
end

% Collect a couple of useful information before starting the actual
% extraction of data and headers
TwixHeader.SiemensVersion       = twix_obj.image.softwareVersion; % Siemens software version (VA,VB,VC,VD,VE?)
TwixHeader.sequenceFileName     = twix_obj.hdr.Config.SequenceFileName; % Full sequence name
TwixHeader.sequenceString       = twix_obj.hdr.Config.SequenceString; % Short sequence name

% Determine the type
% Read information from .image part of the TWIX object
TwixHeader.sqzSize              = twix_obj.image.sqzSize; % dimensions (data points, averages, number of coils, dynamics (ON and OFF))
TwixHeader.sqzDims              = twix_obj.image.sqzDims; % variable names for dimensions
TwixHeader.pointsBeforeEcho     = twix_obj.image.freeParam(1);
TwixData                        = squeeze(twix_obj.image()); % FID data, remove singleton dimensions

% Read information from .hdr part of the TWIX object
TwixHeader.readoutOSFactor      = twix_obj.hdr.Config.ReadoutOSFactor; % Data are oversampled by this factor compared to exam card setting
TwixHeader.removeOS             = twix_obj.hdr.Config.RemoveOversampling; % Is the oversampling removed in the RDA files?
TwixHeader.TR                   = twix_obj.hdr.Config.TR * 1e-3; % TR [ms]
TwixHeader.vectorSize           = twix_obj.hdr.Config.VectorSize; % Data points specified on exam card
TwixHeader.VoI_InPlaneRot       = twix_obj.hdr.Config.VoI_InPlaneRotAngle; % Voxel rotation in plane
TwixHeader.VoI_RoFOV            = twix_obj.hdr.Config.VoI_RoFOV; % Voxel size in readout direction
TwixHeader.VoI_PeFOV            = twix_obj.hdr.Config.VoI_PeFOV; % Voxel size in phase encoding direction
TwixHeader.VoIThickness         = twix_obj.hdr.Config.VoI_SliceThickness; % Voxel size in slice selection direction
TwixHeader.NormCor              = twix_obj.hdr.Config.VoI_Normal_Cor; % Coronal component of normal vector of voxel
TwixHeader.NormSag              = twix_obj.hdr.Config.VoI_Normal_Sag; % Sagittal component of normal vector of voxel
TwixHeader.NormTra              = twix_obj.hdr.Config.VoI_Normal_Tra; % Transversal component of normal vector of voxel
TwixHeader.PosCor               = twix_obj.hdr.Config.VoI_Position_Cor; % Coronal coordinate of voxel
TwixHeader.PosSag               = twix_obj.hdr.Config.VoI_Position_Sag; % Sagittal coordinate of voxel
TwixHeader.PosTra               = twix_obj.hdr.Config.VoI_Position_Tra; % Transversal coordinate of voxel
TwixHeader.SiemensSoftwareVersion  = twix_obj.hdr.Dicom.SoftwareVersions; % Full software version
TwixHeader.Bo                   = twix_obj.hdr.Dicom.flMagneticFieldStrength; % Nominal B0 [T]
TwixHeader.tx_freq              = twix_obj.hdr.Dicom.lFrequency * 1e-6; % Transmitter frequency [MHz]
if iscell(twix_obj.hdr.MeasYaps.alTE)
    TwixHeader.TE               = twix_obj.hdr.MeasYaps.alTE{1} * 1e-3; % TE [ms]
elseif isstruct(twix_obj.hdr.MeasYaps.alTE)
    TwixHeader.TE               = twix_obj.hdr.MeasYaps.alTE(1) * 1e-3; % TE [ms]
end
if iscell(twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime)
    TwixHeader.dwellTime        = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9; % dwell time [s]
elseif isstruct(twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime)
    TwixHeader.dwellTime        = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime(1) * 1e-9; % dwell time [s]
end

% these may only be extractable from a few MEGA-PRESS versions
% editing pulse parameters
if isfield(twix_obj.hdr.MeasYaps, 'sWipMemBlock')
    if isfield(twix_obj.hdr.MeasYaps.sWipMemBlock, 'adFree')
        param = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree;
        param = param(~cellfun('isempty',param));
        TwixHeader.editRF.freq = [param{1}, param{3}+(param{3}-param{1})];
        TwixHeader.editRF.centerFreq = param{3};
        TwixHeader.editRF.bw = param{2};
    end
elseif isfield(twix_obj.hdr.MeasYaps, 'sWiPMemBlock')
    if isfield(twix_obj.hdr.MeasYaps.sWiPMemBlock, 'adFree')
        param = twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree;
        param = param(~cellfun('isempty',param));
        TwixHeader.editRF.freq = [param{1}, param{3}+(param{3}-param{1})];
        TwixHeader.editRF.centerFreq = param{3};
        TwixHeader.editRF.bw = param{2};
    end
end
% delta frequency (center of slice selection)
if isfield(twix_obj.hdr.MeasYaps.sSpecPara, 'dDeltaFrequency')
    TwixHeader.deltaFreq = twix_obj.hdr.MeasYaps.sSpecPara.dDeltaFrequency;
else
    TwixHeader.deltaFreq = 0;
end

% Determine the origin of the sequence
if strfind(TwixHeader.sequenceFileName,'svs_edit')
    TwixHeader.seqtype = 'MEGAPRESS';
    TwixHeader.seqorig = 'WIP'; % Siemens WIP
elseif strfind(TwixHeader.sequenceFileName,'jn_')
    TwixHeader.seqtype = 'MEGAPRESS';
    TwixHeader.seqorig = 'JN'; % Jamie Near's sequence
elseif strfind(TwixHeader.sequenceFileName,'eja_svs_mpress')
    TwixHeader.seqtype = 'MEGAPRESS';
    TwixHeader.seqorig = 'CMRR'; % Minnesota sequence
elseif strfind(TwixHeader.sequenceFileName,'svs_se')
    TwixHeader.seqtype = 'PRESS'; % In case PRESS is used as water reference
else
    TwixHeader.seqorig = TwixHeader.sequenceString;
    error(['Unknown sequence: ' TwixHeader.seqorig '. Please consult the Gannet team for support.'])
end

% Now reorder the FID data array according to software version and sequence 
% origin and sequence type.

if strcmp(TwixHeader.seqtype,'PRESS')
    % For PRESS data, the first dimension of the 4D data array contains the
    % time-domain FID datapoints. The second dimension contains the number
    % of the coils. The third dimension contains the number of averages.
    % The fourth dimension is not well understood, but the second row of
    % this dimension contains all averages, while the first one is empty
    % for all averages but the first one.
    dims.points = 1;
    dims.coils = 2;
    dims.averages = 3;
    dims.dyn = 4;
    TwixData = permute(TwixData,[dims.coils dims.points dims.dyn dims.averages]);
    TwixData = reshape(TwixData,[size(TwixData,1) size(TwixData,2) size(TwixData,3)*size(TwixData,4)]);
elseif strcmp(TwixHeader.seqtype,'MEGAPRESS')    
    % For all known MEGA-PRESS implementations, the first dimension of the 4D
    % data array contains the time-domain FID datapoints.
    dims.points = 1;
    % For all known MEGA-PRESS implementations, the second dimension of the 4D
    % data array contains the the number of the coils.
    dims.coils = 2;
    % It is more difficult for the dimension that contains the averages.
    if strcmp(TwixHeader.SiemensVersion,'vb')
        dims.averages=find(strcmp(TwixHeader.sqzDims,'Set'));
    else
        if strcmp(TwixHeader.seqorig,'CMRR')
            dims.averages=find(strcmp(TwixHeader.sqzDims,'Set'));
        else
            dims.averages=find(strcmp(TwixHeader.sqzDims,'Ave'));
        end
    end
    % It is more difficult for the dimension that contains the dynamics.
    if strcmp(TwixHeader.SiemensVersion,'vb')
        if strcmp(TwixHeader.seqorig,'JN')
            dims.dyn=find(strcmp(TwixHeader.sqzDims,'Ida'));
        else
            dims.dyn=find(strcmp(TwixHeader.sqzDims,'Eco'));
        end
    else
        if strcmp(TwixHeader.seqorig,'CMRR')
            dims.dyn=find(strcmp(TwixHeader.sqzDims,'Eco'));
        elseif strcmp(TwixHeader.seqorig,'JN')
            dims.dyn=find(strcmp(TwixHeader.sqzDims,'Set'));
        else
            dims.dyn=find(strcmp(TwixHeader.sqzDims,'Ide'));
        end
    end
    
    TwixData = permute(TwixData,[dims.coils dims.points dims.dyn dims.averages]);
    TwixData = reshape(TwixData,[size(TwixData,1) size(TwixData,2) size(TwixData,3)*size(TwixData,4)]);
end

end