function MRS_struct = GERead(MRS_struct, fname)
% RTN edits to accommodate Noeske version RAEE 141007
% 160916: MM & RTN edits to accommodate different encoding schemes

ii = MRS_struct.ii;

fid = fopen(fname, 'r', 'ieee-be');
if fid == -1
    tmp = ['Unable to locate Pfile ' fname];
    disp(tmp);
    return
end
% Return error message if unable to read file type.
% Determine size of Pfile header based on rdbm_rev number
fseek(fid, 0, 'bof');
f_hdr_value = fread(fid, 1, 'real*4');
rdbm_rev_num = f_hdr_value(1);
if rdbm_rev_num == 7.0
    pfile_header_size = 39984; % LX
elseif rdbm_rev_num == 8.0
    pfile_header_size = 60464; % Cardiac / MGD
elseif rdbm_rev_num > 5.0 && rdbm_rev_num < 6.0
    pfile_header_size = 39940; % Signa 5.5
else
    % In 11.0 and later the header and data are stored as little-endian
    fclose(fid);
    fid = fopen(fname, 'r', 'ieee-le');
    fseek(fid, 0, 'bof');
    f_hdr_value = fread(fid, 1, 'real*4');
    if f_hdr_value == 9.0 % 11.0 product release
        pfile_header_size = 61464;
    elseif f_hdr_value == 11.0 % 12.0 product release
        pfile_header_size = 66072;
    elseif f_hdr_value > 11.0 && f_hdr_value < 100.0 % 14.0 and later
        fseek(fid, 1468, 'bof');
        pfile_header_size = fread(fid, 1, 'integer*4');
    else
        sprintf('Invalid Pfile header revision: %f', f_hdr_value);
        return
    end
end

MRS_struct.p.GE.rdbm_rev_num = f_hdr_value(1); % MM (170118)
chkRev = {'14.3', '16', '24'}; % GERead mods tested with these revisions only
if ~any(strcmp(num2str(MRS_struct.p.GE.rdbm_rev_num), chkRev))
    warning('GERead not fully functional with header revision number %g!', MRS_struct.p.GE.rdbm_rev_num);
end
% Read header information
fseek(fid, 0, 'bof');
hdr_value = fread(fid, 102, 'integer*2');
% RTN - read rhuser
fseek(fid, 0, 'bof');
f_hdr_value = fread(fid, 74, 'real*4');
% RTN (170118): Find center frequency
fseek(fid, 0, 'bof');
i_hdr_value = fread(fid, 102+9, 'integer*4');
MRS_struct.p.LarmorFreq(ii) = i_hdr_value(102+5)/1e7;
MRS_struct.p.sw(ii) = f_hdr_value(55); % MM (160916)

nechoes = hdr_value(36);
MRS_struct.p.GE.nechoes = nechoes;
% RTN - number of phase cycles
nex = hdr_value(37);
MRS_struct.p.GE.NEX = nex;
nframes = hdr_value(38);
point_size = hdr_value(42);
MRS_struct.p.npoints(ii) = hdr_value(52);
MRS_struct.p.nrows(ii) = hdr_value(53);
start_recv = hdr_value(101);
stop_recv = hdr_value(102);
nreceivers = (stop_recv - start_recv) + 1;

% MM (170118): Find TE/TR
fseek(fid, 1468, 'bof');
p_hdr_value = fread(fid, 12, 'integer*4'); % byte offsets to start of sub-header structures
fseek(fid, p_hdr_value(10), 'bof'); % set position to start of rdb_hdr_image
t_hdr_value = fread(fid, p_hdr_value(1)-p_hdr_value(10), 'integer*4');
if strcmp(num2str(MRS_struct.p.GE.rdbm_rev_num), '14.3')
    MRS_struct.p.TE(ii) = t_hdr_value(181)/1e3;
    MRS_struct.p.TR(ii) = t_hdr_value(179)/1e3;
elseif MRS_struct.p.GE.rdbm_rev_num == 16
    MRS_struct.p.TE(ii) = t_hdr_value(193)/1e3;
    MRS_struct.p.TR(ii) = t_hdr_value(191)/1e3;
elseif MRS_struct.p.GE.rdbm_rev_num == 24
    MRS_struct.p.TE(ii) = t_hdr_value(267)/1e3;
    MRS_struct.p.TR(ii) = t_hdr_value(265)/1e3;
end

% MM (170127): Find voxel dimensions and edit pulse parameters
fseek(fid, p_hdr_value(8), 'bof'); % set position to start of rdb_hdr_exam
o_hdr_value = fread(fid, p_hdr_value(9)-p_hdr_value(8), 'real*4');
if strcmp(num2str(MRS_struct.p.GE.rdbm_rev_num), '14.3')
    MRS_struct.p.voxdim(ii,:) = o_hdr_value(810:812)';
    MRS_struct.p.GE.editRF.waveform(ii) = o_hdr_value(821);
    MRS_struct.p.GE.editRF.freq_Hz(ii,:) = o_hdr_value(822:823)';
    MRS_struct.p.GE.editRF.freq_ppm(ii,:) = (MRS_struct.p.GE.editRF.freq_Hz(ii,:) / MRS_struct.p.LarmorFreq(ii)) + 4.68;
    MRS_struct.p.GE.editRF.dur(ii) = o_hdr_value(824)/1e3;
elseif MRS_struct.p.GE.rdbm_rev_num == 16
    MRS_struct.p.voxdim(ii,:) = o_hdr_value(822:824)';
    MRS_struct.p.GE.editRF.waveform(ii) = o_hdr_value(833);
    MRS_struct.p.GE.editRF.freq_Hz(ii,:) = o_hdr_value(834:835)';
    MRS_struct.p.GE.editRF.freq_ppm(ii,:) = (MRS_struct.p.GE.editRF.freq_Hz(ii,:) / MRS_struct.p.LarmorFreq(ii)) + 4.68;
    MRS_struct.p.GE.editRF.dur(ii) = o_hdr_value(836)/1e3;
elseif MRS_struct.p.GE.rdbm_rev_num == 24
    MRS_struct.p.voxdim(ii,:) = o_hdr_value(1228:1230)';
    MRS_struct.p.GE.editRF.waveform(ii) = o_hdr_value(1239);
    MRS_struct.p.GE.editRF.freq_Hz(ii,:) = o_hdr_value(1240:1241)';
    MRS_struct.p.GE.editRF.freq_ppm(ii,:) = (MRS_struct.p.GE.editRF.freq_Hz(ii,:) / MRS_struct.p.LarmorFreq(ii)) + 4.68;
    MRS_struct.p.GE.editRF.dur(ii) = o_hdr_value(1242)/1e3;
end

% Spectro prescan pfiles
if MRS_struct.p.npoints(ii) == 1 && MRS_struct.p.nrows(ii) == 1
    MRS_struct.p.npoints(ii) = 2048;
end

% Compute size (in bytes) of each frame, echo and slice
data_elements = MRS_struct.p.npoints(ii)*2;
frame_size = data_elements*point_size;
my_frame = 1;

% Start to read data into Eightchannel structure.
totalframes = (MRS_struct.p.nrows(ii)-my_frame+1)*nechoes; % RTN nechoes mulitply
MRS_struct.p.nrows(ii) = totalframes;
data_elements2 = data_elements * totalframes * nreceivers;

% Compute offset in bytes to start of frame.
file_offset = pfile_header_size + ((my_frame-1)*frame_size);

fseek(fid, file_offset, 'bof');
% Read data: point_size = 2 means 16 bit data, point_size = 4 means EDR
if point_size == 2
    raw_data = fread(fid, data_elements2, 'integer*2');
else
    raw_data = fread(fid, data_elements2, 'integer*4');
end
fclose(fid);

% 110303 CJE
% Calculate Navg from nframes, 8 water frames, 2 phase cycles
% Needs to be specific to single experiment - for frame rejection
% RTN edits to accommodate Noeske version raee 20141007
% MM (160916): Incorporating more edits from RTN to handle dual-echo data
%              acquired with one of four possible encoding schemes:
%              NEX=2/noadd=0, NEX=2/noadd=1, NEX=8/noadd=0, NEX=8/noadd=1
% MM (171120): RTN edits to accomodate HERMES aquisitions; better looping
%              over phase cycles
if nechoes == 1
    MRS_struct.p.Navg(ii) = (nframes-8)*2;
    MRS_struct.p.Nwateravg = 8;
    ShapeData = reshape(raw_data,[2 MRS_struct.p.npoints(ii) totalframes nreceivers]);
    WaterData = ShapeData(:,:,2:9,:);
    FullData = ShapeData(:,:,10:end,:);
    
    totalframes = totalframes-9;
    MRS_struct.p.nrows(ii)=totalframes;
    
    Frames_for_Water = 8;
else
    dataframes = f_hdr_value(59)/nex;
    refframes = f_hdr_value(74);
    
    MRS_struct.p.Navg(ii) = dataframes * nex * nechoes; % RTN 2017
    
    if (dataframes+refframes) ~= nframes
        mult = nex/2; % RTN 2016
        multw = nex; % RTN 2016
        %mult = 1; % RTN 2017
        %multw = 1; % RTN 2017
        MRS_struct.p.GE.noadd = 1;
        dataframes = dataframes * nex;
        refframes = nframes - dataframes; % refframes*nex; 2015
    else
        mult = nex/2; % RTN 2016
        multw = 1; % RTN 2016
        %mult = 1; % RTN 2017
        %multw = 1/nex; % RTN 2017
        MRS_struct.p.GE.noadd = 0;
    end
    
    MRS_struct.p.Nwateravg(ii) = refframes * nechoes; % RTN 2017
    
    if totalframes ~= (dataframes+refframes+1)*nechoes % RTN 2017
        error('# of totalframes not same as (dataframes+refframes+1)*2');
    end
    ShapeData = reshape(raw_data,[2 MRS_struct.p.npoints(ii) totalframes nreceivers]);
    WaterData = zeros([2 MRS_struct.p.npoints(ii) refframes*nechoes nreceivers]); % RTN 2017
    for loop = 1:refframes
        for echoloop = 1:nechoes
            WaterData(:,:,nechoes*(loop-1)+echoloop,:) = (-1)^(MRS_struct.p.GE.noadd*(loop-1)) * ...
                ShapeData(:,:,1+(totalframes/nechoes)*(echoloop-1)+loop,:) * multw; % RTN 2017
        end
    end
    FullData = zeros([2 MRS_struct.p.npoints(ii) dataframes*2 nreceivers]);
    for loop = 1:dataframes
        for echoloop = 1:nechoes
            FullData(:,:,nechoes*(loop-1)+echoloop,:)=(-1)^(MRS_struct.p.GE.noadd*(loop-1)) * ...
                ShapeData(:,:,1+refframes+(totalframes/nechoes)*(echoloop-1)+loop,:) * mult; % RTN 2017
        end
    end
    totalframes = totalframes - (refframes+1) * nechoes; % RTN 2017
    MRS_struct.p.nrows(ii) = totalframes;
    Frames_for_Water = refframes * nechoes; % RTN 2017
end

FullData = FullData.*repmat([1;1i], [1 MRS_struct.p.npoints(ii) totalframes nreceivers]);
FullData = squeeze(sum(FullData,1));
FullData = permute(FullData,[3 1 2]);
WaterData = WaterData.*repmat([1;1i], [1 MRS_struct.p.npoints(ii) Frames_for_Water nreceivers]);
WaterData = squeeze(sum(WaterData,1));
WaterData = permute(WaterData,[3 1 2]);

% MM (170505)
firstpoint_water = conj(WaterData(:,1,:));
channels_scale = squeeze(sqrt(sum(firstpoint_water .* conj(firstpoint_water),1)));
channels_scale = repmat(channels_scale, [1 nreceivers MRS_struct.p.npoints(ii)]);
channels_scale = permute(channels_scale, [2 3 1]);
firstpoint_water = repmat(firstpoint_water, [1 MRS_struct.p.npoints(ii) 1])./channels_scale;

WaterData = WaterData .* firstpoint_water;
WaterData = squeeze(sum(WaterData,1));
MRS_struct.fids.data_water = WaterData;

% Use first point of water data to phase water-suppressed data
firstpoint = mean(firstpoint_water,3);
firstpoint = repmat(firstpoint, [1 1 size(FullData,3)]);

FullData = FullData .* firstpoint;
FullData = squeeze(sum(FullData,1));
MRS_struct.fids.data = FullData;

rescale = 1/1e11; % necessary for GE data or numbers blow up
MRS_struct.fids.data = MRS_struct.fids.data * rescale;
MRS_struct.fids.data_water = MRS_struct.fids.data_water * rescale;

end



