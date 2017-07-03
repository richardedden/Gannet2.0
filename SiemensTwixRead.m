function MRS_struct = SiemensTwixRead(MRS_struct, fname, fname_water)
%131216 Since twix data is not combined, use same code from GERead to bring
%in Siemens twix data

%TWIX data can either come from Jamie Near's sequence or the
%Siemens WIP.  Jamie's has NSet as 2 and the WIP has NEco.

ii = MRS_struct.ii;

%This handles the GABA data - it is needed whatever..
%Use mapVBVD to pull in data.
twix_obj=mapVBVD(fname);
if any(MRS_struct.p.Siemens_type == [4 5 6 7])
    twix_obj=twix_obj{2};
end
pointsBeforeEcho = twix_obj.image.freeParam(1);

%This code included by kind permission of Jamie Near.
%Pull in some header information not accessed by mapVBVD

% Get spectral width and dwelltime
% Check if twix_obj header is available
if isfield(twix_obj,'hdr')
    dwelltime=twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1} * 1e-9;
    spectralwidth=1/dwelltime;
else
    fid=fopen(fname);
    line=fgets(fid);
    index=strfind(line,'sRXSPEC.alDwellTime[0]');
    equals_index=strfind(line,'= ');
    while isempty(index) || isempty(equals_index)
        line=fgets(fid);
        index=strfind(line,'sRXSPEC.alDwellTime[0]');
        equals_index=strfind(line,'= ');
    end
    dwelltime=line(equals_index+1:end);
    dwelltime=str2double(dwelltime)*1e-9;
    spectralwidth=1/dwelltime;
    fclose(fid);
end
%End of Jamie Near's code

%Calculate some parameters:
MRS_struct.p.sw = spectralwidth;
MRS_struct.p.LarmorFreq = twix_obj.hdr.Config.Frequency/1e6; % MM (170104)
MRS_struct.p.nrows = twix_obj.image.NAcq;

% MM (170127)
MRS_struct.p.TR(ii) = twix_obj.hdr.Phoenix.alTR{1}/1e3;
MRS_struct.p.TE(ii) = twix_obj.hdr.Phoenix.alTE{1}/1e3;
MRS_struct.p.npoints = twix_obj.image.NCol;
if isfield(twix_obj.hdr.Meas, 'VoI_RoFOV')
    MRS_struct.p.voxdim(ii,:) = [twix_obj.hdr.Meas.VoI_RoFOV ...
        twix_obj.hdr.Meas.VoI_PeFOV ...
        twix_obj.hdr.Meas.VoiThickness];
else
    MRS_struct.p.voxdim(ii,:) = [twix_obj.hdr.MeasYaps.sSpecPara.sVoI.dThickness ...
        twix_obj.hdr.MeasYaps.sSpecPara.sVoI.dPhaseFOV ...
        twix_obj.hdr.MeasYaps.sSpecPara.sVoI.dReadoutFOV];
end
if isfield(twix_obj.hdr.MeasYaps, 'sWipMemBlock')
    param = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree;
else
    param = twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree;
end
param = param(~cellfun('isempty',param));
MRS_struct.p.Siemens.editRF.freq(ii,:) = [param{1}, param{3}+(param{3}-param{1})];
MRS_struct.p.Siemens.editRF.centerFreq(ii) = param{3};
MRS_struct.p.Siemens.editRF.bw(ii) = param{2};
if isfield(twix_obj.hdr.MeasYaps.sSpecPara, 'dDeltaFrequency')
    MRS_struct.p.Siemens.deltaFreq.metab(ii) = twix_obj.hdr.MeasYaps.sSpecPara.dDeltaFrequency;
else
    MRS_struct.p.Siemens.deltaFreq.metab(ii) = 0;
end

switch MRS_struct.p.Siemens_type
    case 1
        %[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NEco twix_obj.image.NSet];
        FullData=permute(reshape(double(twix_obj.image()),[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NEco twix_obj.image.NSet]),[2 1 3 4]);
        %Undo Plus-minus
        %FullData(:,:,2,:)=-FullData(:,:,2,:);
        FullData=reshape(FullData,[twix_obj.image.NCha twix_obj.image.NCol twix_obj.image.NSet*twix_obj.image.NEco]);
    case 2
        FullData=permute(reshape(double(twix_obj.image()),[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NSet twix_obj.image.NIda]),[2 1 4 3]);
        %Undo Plus-minus
        FullData(:,:,2,:)=-FullData(:,:,2,:);
        FullData=reshape(FullData,[twix_obj.image.NCha twix_obj.image.NCol twix_obj.image.NSet*twix_obj.image.NIda]);
    case 3
        %[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NAve twix_obj.image.NIde];
        FullData=permute(reshape(double(twix_obj.image()),[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NAve twix_obj.image.NIde]),[2 1 4 3]);
        %Undo Plus-minus
        %FullData(:,:,2,:)=-FullData(:,:,2,:);
        %size(FullData)
        FullData=reshape(FullData,[twix_obj.image.NCha twix_obj.image.NCol twix_obj.image.NAve*twix_obj.image.NIde]);
    case 4
        %size(twix_obj.image())
        %[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NAve twix_obj.image.NIde];
        %[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NEco twix_obj.image.NSet];
        FullData=permute(reshape(double(twix_obj.image()),[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NAve twix_obj.image.NIde]),[2 1 4 3]);
        %Undo Plus-minus
        %FullData(:,:,2,:)=-FullData(:,:,2,:);
        %size(FullData)
        FullData=reshape(FullData,[twix_obj.image.NCha twix_obj.image.NCol twix_obj.image.NAve*twix_obj.image.NIde]);
    case 5
        %size(twix_obj.image())
        %[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NAve twix_obj.image.NIde];
        %[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NEco twix_obj.image.NSet];
        FullData=permute(reshape(double(twix_obj.image()),[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NAve twix_obj.image.NIde]),[2 1 4 3]);
        size(FullData)
        %Undo Plus-minus
        %FullData(:,:,2,:)=-FullData(:,:,2,:);
        %size(FullData)
        FullData=reshape(FullData,[twix_obj.image.NCha twix_obj.image.NCol twix_obj.image.NAve*twix_obj.image.NIde]);
    case 6
        FullData=permute(reshape(double(twix_obj.image()),[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NAve twix_obj.image.NSet]),[2 1 4 3]);
        FullData=reshape(FullData,[twix_obj.image.NCha twix_obj.image.NCol twix_obj.image.NAve*twix_obj.image.NSet]);
    case 7
        %[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NEco twix_obj.image.NSet];
        FullData=permute(reshape(double(twix_obj.image()),[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NAve twix_obj.image.NSet]),[2 1 4 3]);
        FullData=reshape(FullData,[twix_obj.image.NCha twix_obj.image.NCol twix_obj.image.NAve*twix_obj.image.NSet]);
end

MRS_struct.p.Navg(ii) = double(twix_obj.image.NAcq);
FullData = FullData(:,pointsBeforeEcho+1:end,:);
MRS_struct.p.npoints = MRS_struct.p.npoints - pointsBeforeEcho; % MM (160914)

% If no water data provided, combine data based upon first point of
% water-suppressed data (average all transients) (MM: 170124)
if nargin < 3
    
    firstpoint = mean(conj(FullData(:,1,:)),3);
    channels_scale = squeeze(sqrt(sum(firstpoint.*conj(firstpoint))));
    firstpoint = repmat(firstpoint, [1 MRS_struct.p.npoints MRS_struct.p.nrows])/channels_scale;
    FullData = FullData .* firstpoint;
    FullData = conj(squeeze(sum(FullData,1)));
    MRS_struct.fids.data = FullData;
    
end

if nargin > 2 % water data provided
    
    %Then we additionally need to pull in the water data.
    twix_obj_water=mapVBVD(fname_water);
    if any(MRS_struct.p.Siemens_type == [4 5 6 7])
        twix_obj_water=twix_obj_water{2};
    end
    pointsBeforeEcho = twix_obj_water.image.freeParam(1); % MM (160914)
    MRS_struct.p.nrows_water = twix_obj_water.image.NAcq;
    MRS_struct.p.npoints_water = twix_obj_water.image.NCol;
    MRS_struct.p.Nwateravg(ii) = MRS_struct.p.nrows_water;
    
    % MM (170127)
    if isfield(twix_obj_water.hdr.MeasYaps.sSpecPara,'dDeltaFrequency')
        MRS_struct.p.Siemens.deltaFreq.water(ii) = twix_obj_water.hdr.MeasYaps.sSpecPara.dDeltaFrequency;
    else
        MRS_struct.p.Siemens.deltaFreq.water(ii) = 0;
    end
    MRS_struct.p.Siemens = reorderstructure(MRS_struct.p.Siemens, 'editRF', 'deltaFreq');
    
    if twix_obj.image.NEco > twix_obj.image.NIda
        
        WaterData = permute(reshape(double(twix_obj_water.image()), [twix_obj_water.image.NCol, twix_obj_water.image.NCha, twix_obj_water.image.NEco, twix_obj_water.image.NSet]), [2 1 3 4]);
        WaterData = reshape(WaterData, [twix_obj_water.image.NCha, twix_obj_water.image.NCol, twix_obj_water.image.NSet*twix_obj_water.image.NEco]);
        
    else
        
        % added below to make work for Skyra, other cases not tested
        % and are copies of what was there originally was there but
        % according to above, in some the undo plus-mins on the
        % WaterData(:,:,2,:) may not be needed
        switch MRS_struct.p.Siemens_type
            case 1 % this was the original code - made it into case 1, 2
                % Copy it into WaterData
                WaterData=permute(reshape(double(twix_obj_water.image()),[twix_obj_water.image.NCol twix_obj_water.image.NCha twix_obj_water.image.NSet twix_obj_water.image.NIda]),[2 1 4 3]);
                %Undo Plus-minus
                WaterData(:,:,2,:)=-WaterData(:,:,2,:);
                WaterData=reshape(WaterData,[twix_obj_water.image.NCha twix_obj_water.image.NCol twix_obj_water.image.NSet*twix_obj_water.image.NIda]);
            case 2 % this was the original code - made it into case 1, 2
                % Copy it into WaterData
                WaterData=permute(reshape(double(twix_obj_water.image()),[twix_obj_water.image.NCol twix_obj_water.image.NCha twix_obj_water.image.NSet twix_obj_water.image.NIda]),[2 1 4 3]);
                %Undo Plus-minus
                WaterData(:,:,2,:)=-WaterData(:,:,2,:);
                WaterData=reshape(WaterData,[twix_obj_water.image.NCha twix_obj_water.image.NCol twix_obj_water.image.NSet*twix_obj_water.image.NIda]);
            case 3 % above in case 3 didn't do the WaterData(:,:,2,:)=-WaterData(:,:,2,:); so commented out here.
                % Copy it into WaterData
                WaterData=permute(reshape(double(twix_obj_water.image()),[twix_obj_water.image.NCol twix_obj_water.image.NCha twix_obj_water.image.NAve twix_obj_water.image.NIde]),[2 1 4 3]);
                WaterData=reshape(WaterData,[twix_obj_water.image.NCha twix_obj_water.image.NCol twix_obj_water.image.NAve*twix_obj_water.image.NIde]);
            case 4 %
                % Copy it into WaterData
                WaterData=permute(reshape(double(twix_obj_water.image()),[twix_obj_water.image.NCol twix_obj_water.image.NCha twix_obj_water.image.NAve twix_obj_water.image.NIde]),[2 1 4 3]);
                WaterData=reshape(WaterData,[twix_obj_water.image.NCha twix_obj_water.image.NCol twix_obj_water.image.NAve*twix_obj_water.image.NIde]);
            case 6 %
                % Copy it into WaterData
                WaterData=permute(reshape(double(twix_obj_water.image()),[twix_obj_water.image.NCol twix_obj_water.image.NCha twix_obj_water.image.NAve twix_obj_water.image.NSet]),[2 1 4 3]);
                %Undo Plus-minus
                WaterData(:,:,2,:)=-WaterData(:,:,2,:);
                WaterData=reshape(WaterData,[twix_obj_water.image.NCha twix_obj_water.image.NCol twix_obj_water.image.NAve*twix_obj_water.image.NSet]);
            case 7
                % Copy it into WaterData
                WaterData=permute(reshape(double(twix_obj_water.image()),[twix_obj_water.image.NCol twix_obj_water.image.NCha twix_obj_water.image.NAve twix_obj_water.image.NSet]),[2 1 4 3]);
                %Undo Plus-minus
                WaterData(:,:,2,:)=-WaterData(:,:,2,:);
                WaterData=reshape(WaterData,[twix_obj_water.image.NCha twix_obj_water.image.NCol twix_obj_water.image.NAve*twix_obj_water.image.NSet]);
        end
        
    end
    
    % MM (160914)
    WaterData = WaterData(:,pointsBeforeEcho+1:end,:);
    MRS_struct.p.npoints_water = MRS_struct.p.npoints_water - pointsBeforeEcho;
    
    % MM (170123)
    firstpoint_water = conj(WaterData(:,1,:));
    channels_scale = squeeze(sqrt(sum(firstpoint_water .* conj(firstpoint_water),1)));
    channels_scale = repmat(channels_scale, [1 twix_obj_water.image.NCha MRS_struct.p.npoints_water]);
    channels_scale = permute(channels_scale, [2 3 1]);
    firstpoint_water = repmat(firstpoint_water, [1 MRS_struct.p.npoints_water 1])./channels_scale;
    WaterData = WaterData .* firstpoint_water;
    WaterData = conj(squeeze(sum(WaterData,1)));
    WaterData = squeeze(mean(WaterData,2));
    MRS_struct.fids.data_water = WaterData;
    
    % Use first point of water data to phase water-suppressed data
    if MRS_struct.p.Siemens_type == 7 % +/- phasing in water FIDs
        firstpoint_water(:,:,2:2:end) = -firstpoint_water(:,:,2:2:end);
        firstpoint = mean(firstpoint_water,3);
    else
        firstpoint = mean(firstpoint_water,3);
    end
    firstpoint = repmat(firstpoint, [1 1 size(FullData,3)]);
    FullData = FullData .* firstpoint;
    FullData = conj(squeeze(sum(FullData,1)));
    MRS_struct.fids.data = FullData;
    
end

end



