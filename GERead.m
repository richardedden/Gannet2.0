function MRS_struct = GERead(MRS_struct, fname)
% 141007: RTN edits to accommodate Noeske version
% 160916: MM & RTN edits to accommodate different encoding schemes
% 180404: RTN edits for more flexible handling of different P-file header
% revisions; added support for rdbm_rev_num 26.002

ii = MRS_struct.ii;

fid = fopen(fname, 'r', 'ieee-be');
fseek(fid, 0, 'bof');
rdbm_rev_num = fread(fid, 1, 'real*4');
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
    rdbm_rev_num = fread(fid, 1, 'real*4');
    if rdbm_rev_num == 9.0 % 11.0 product release
        pfile_header_size = 61464;
    elseif rdbm_rev_num == 11.0 % 12.0 product release
        pfile_header_size = 66072;
    end
end

MRS_struct.p.GE.rdbm_rev_num = rdbm_rev_num;
chkRev = {'14.3', '16', '24', '26.002'}; % GERead mods tested with these revisions only
if ~any(strcmp(num2str(rdbm_rev_num), chkRev))
    warning('GERead not fully functional with header revision number %g!', rdbm_rev_num);
end

% RTN 2018
% Added flexible P-file revision support
% Values are read from rdb_hdr and image sub-headers
% Position can be found in rdbm.h (RDB_HEADER_REC) and imagedb.h (MRIMAGEDATATYPE)

% RTN 2018
% unsigned int rdb_hdr_ps_mps_freq
% float rdb_hdr_user0
% float rdb_hdr_user4
% float rdb_hdr_user19
% short rdb_hdr_nechoes
% short rdb_hdr_navs
% short rdb_hdr_nframes
% short rdb_hdr_point_size
% unsigned short rdb_hdr_da_xres
% short rdb_hdr_da_yres
% short rdb_hdr_dab[0].start_rcv
% short rdb_hdr_dab[0].stop_rcv
% int rdb_hdr_off_image
% int rdb_hdr_off_data
%
% image sub-header
% int te
% int tr
% float user8-10    voxel dimensions
% float user19      rf waveform
% float user20-21   offset frequencies
% float user22      pulse width (-1 default)

switch num2str(rdbm_rev_num)
    
    case '14.3'
        
        % int
        rdb_hdr_off_image   = 377;
        rdb_hdr_off_data    = 368;
        rdb_hdr_ps_mps_freq = 107;
        
        % float
        rdb_hdr_user0  = 55;
        rdb_hdr_user4  = 59;
        rdb_hdr_user19 = 74;
        
        % short
        rdb_hdr_nechoes       = 36;
        rdb_hdr_navs          = 37;
        rdb_hdr_nframes       = 38;
        rdb_hdr_point_size    = 42;
        rdb_hdr_da_xres       = 52;
        rdb_hdr_da_yres       = 53;
        rdb_hdr_dab_start_rcv = 101;
        rdb_hdr_dab_stop_rcv  = 102;
        
        % int
        image_te = 181;
        image_tr = 179;
        
        % float
        image_user8  = 38;
        image_user19 = 49;
        image_user20 = 50;
        image_user22 = 52;
        
    case '16'
        
        % int
        rdb_hdr_off_image   = 377;
        rdb_hdr_off_data    = 368;
        rdb_hdr_ps_mps_freq = 107;
        
        % float
        rdb_hdr_user0  = 55;
        rdb_hdr_user4  = 59;
        rdb_hdr_user19 = 74;
        
        % short
        rdb_hdr_nechoes       = 36;
        rdb_hdr_navs          = 37;
        rdb_hdr_nframes       = 38;
        rdb_hdr_point_size    = 42;
        rdb_hdr_da_xres       = 52;
        rdb_hdr_da_yres       = 53;
        rdb_hdr_dab_start_rcv = 101;
        rdb_hdr_dab_stop_rcv  = 102;
        
        % int
        image_te = 193;
        image_tr = 191;
        
        % float
        image_user8  = 50;
        image_user19 = 61;
        image_user20 = 62;
        image_user22 = 64;
        
    case '24'
        
        % int
        rdb_hdr_off_image   = 377;
        rdb_hdr_off_data    = 368;
        rdb_hdr_ps_mps_freq = 107;
        
        % float
        rdb_hdr_user0  = 55;
        rdb_hdr_user4  = 59;
        rdb_hdr_user19 = 74;
        
        % short
        rdb_hdr_nechoes       = 36;
        rdb_hdr_navs          = 37;
        rdb_hdr_nframes       = 38;
        rdb_hdr_point_size    = 42;
        rdb_hdr_da_xres       = 52;
        rdb_hdr_da_yres       = 53;
        rdb_hdr_dab_start_rcv = 101;
        rdb_hdr_dab_stop_rcv  = 102;
        
        % int
        image_te = 267;
        image_tr = 265;
        
        % float
        image_user8  = 98;
        image_user19 = 109;
        image_user20 = 110;
        image_user22 = 112;
        
    case '26.002'
        
        % int
        rdb_hdr_off_image   = 11;
        rdb_hdr_off_data    = 2;
        rdb_hdr_ps_mps_freq = 123;
        
        % float
        rdb_hdr_user0  = 71;
        rdb_hdr_user4  = 75;
        rdb_hdr_user19 = 90;
        
        % short
        rdb_hdr_nechoes       = 74;
        rdb_hdr_navs          = 75;
        rdb_hdr_nframes       = 76;
        rdb_hdr_point_size    = 80;
        rdb_hdr_da_xres       = 90;
        rdb_hdr_da_yres       = 91;
        rdb_hdr_dab_start_rcv = 133;
        rdb_hdr_dab_stop_rcv  = 134;
        
        % int
        image_te = 267;
        image_tr = 265;
        
        % float
        image_user8  = 98;
        image_user19 = 109;
        image_user20 = 110;
        image_user22 = 112;
        
end

% Read rdb header as short, int and float
fseek(fid, 0, 'bof');
hdr_value = fread(fid, rdb_hdr_dab_stop_rcv, 'integer*2');
fseek(fid, 0, 'bof');
f_hdr_value = fread(fid, rdb_hdr_user19, 'real*4');
fseek(fid, 0, 'bof');
i_hdr_value = fread(fid, max(rdb_hdr_off_image, rdb_hdr_ps_mps_freq), 'integer*4');

if rdbm_rev_num > 11.0
    pfile_header_size = i_hdr_value(rdb_hdr_off_data);
end

MRS_struct.p.LarmorFreq(ii) = i_hdr_value(rdb_hdr_ps_mps_freq)/1e7;
MRS_struct.p.sw(ii) = f_hdr_value(rdb_hdr_user0);

nechoes = hdr_value(rdb_hdr_nechoes);
MRS_struct.p.GE.nechoes = nechoes;
nex = hdr_value(rdb_hdr_navs);
MRS_struct.p.GE.NEX = nex;
nframes = hdr_value(rdb_hdr_nframes);
point_size = hdr_value(rdb_hdr_point_size);
MRS_struct.p.npoints(ii) = hdr_value(rdb_hdr_da_xres);
MRS_struct.p.nrows(ii) = hdr_value(rdb_hdr_da_yres);

start_recv = hdr_value(rdb_hdr_dab_start_rcv);
stop_recv = hdr_value(rdb_hdr_dab_stop_rcv);
nreceivers = (stop_recv - start_recv) + 1;

% RTN 2018
dataframes = f_hdr_value(rdb_hdr_user4)/nex;
refframes = f_hdr_value(rdb_hdr_user19);

% Read image header as int and float
% MM (170118): Find TE/TR
fseek(fid, i_hdr_value(rdb_hdr_off_image), 'bof');
t_hdr_value = fread(fid, image_te, 'integer*4');
fseek(fid, i_hdr_value(rdb_hdr_off_image), 'bof');
o_hdr_value = fread(fid, image_user22, 'real*4');
MRS_struct.p.TE(ii) = t_hdr_value(image_te)/1e3;
MRS_struct.p.TR(ii) = t_hdr_value(image_tr)/1e3;

% MM (170127): Find voxel dimensions and edit pulse parameters
MRS_struct.p.voxdim(ii,:) = o_hdr_value(image_user8:image_user8+2)';
MRS_struct.p.GE.editRF.waveform(ii) = o_hdr_value(image_user19);
MRS_struct.p.GE.editRF.freq_Hz(ii,:) = o_hdr_value(image_user20:image_user20+1)';
MRS_struct.p.GE.editRF.freq_ppm(ii,:) = (MRS_struct.p.GE.editRF.freq_Hz(ii,:) / MRS_struct.p.LarmorFreq(ii)) + 4.68;
MRS_struct.p.GE.editRF.dur(ii) = o_hdr_value(image_user22)/1e3;
% RTN 2018: check for default value (-1) of pulse length
if MRS_struct.p.GE.editRF.dur(ii) <= 0
    MRS_struct.p.GE.editRF.dur(ii) = 16;
end

% Spectro prescan pfiles
if MRS_struct.p.npoints(ii) == 1 && MRS_struct.p.nrows(ii) == 1
    MRS_struct.p.npoints(ii) = 2048;
end

% Compute size (in bytes) of data
data_elements = MRS_struct.p.npoints(ii) * 2;
totalframes = MRS_struct.p.nrows(ii) * nechoes; % RTN nechoes mulitply
MRS_struct.p.nrows(ii) = totalframes;
data_elements = data_elements * totalframes * nreceivers;

fseek(fid, pfile_header_size, 'bof');
% Read data: point_size = 2 means 16-bit data, point_size = 4 means EDR
if point_size == 2
    raw_data = fread(fid, data_elements, 'integer*2');
else
    raw_data = fread(fid, data_elements, 'integer*4');
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
    
    MRS_struct.p.Navg(ii) = (nframes - 8) * 2;
    MRS_struct.p.Nwateravg = 8;
    ShapeData = reshape(raw_data, [2 MRS_struct.p.npoints(ii) totalframes nreceivers]);
    WaterData = ShapeData(:,:,2:9,:);
    FullData = ShapeData(:,:,10:end,:);
    
    totalframes = totalframes - 9;
    MRS_struct.p.nrows(ii) = totalframes;
    waterframes = 8;
    
else
    
    MRS_struct.p.Navg(ii) = dataframes * nex * nechoes; % RTN 2017
    
    if (dataframes + refframes) ~= nframes
        mult = nex/2; % RTN 2016
        multw = nex; % RTN 2016
        %mult = 1; % RTN 2017
        %multw = 1; % RTN 2017
        MRS_struct.p.GE.noadd = 1;
        dataframes = dataframes * nex;
        refframes = nframes - dataframes;
    else
        mult = nex/2; % RTN 2016
        multw = 1; % RTN 2016
        %mult = 1; % RTN 2017
        %multw = 1/nex; % RTN 2017
        MRS_struct.p.GE.noadd = 0;
    end
    
    MRS_struct.p.Nwateravg(ii) = refframes * nechoes; % RTN 2017
    
    if totalframes ~= (dataframes + refframes + 1) * nechoes % RTN 2017
        error('# of totalframes not same as (dataframes + refframes + 1) * nechoes');
    end
    ShapeData = reshape(raw_data, [2 MRS_struct.p.npoints(ii) totalframes nreceivers]);
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
            FullData(:,:,nechoes*(loop-1)+echoloop,:) = (-1)^(MRS_struct.p.GE.noadd*(loop-1)) * ...
                ShapeData(:,:,1+refframes+(totalframes/nechoes)*(echoloop-1)+loop,:) * mult; % RTN 2017
        end
    end
    totalframes = totalframes - (refframes + 1) * nechoes; % RTN 2017
    MRS_struct.p.nrows(ii) = totalframes;
    waterframes = refframes * nechoes; % RTN 2017
end

FullData = FullData .* repmat([1; 1i], [1 MRS_struct.p.npoints(ii) totalframes nreceivers]);
FullData = squeeze(sum(FullData,1));
FullData = permute(FullData, [3 1 2]);
WaterData = WaterData .* repmat([1; 1i], [1 MRS_struct.p.npoints(ii) waterframes nreceivers]);
WaterData = squeeze(sum(WaterData,1));
WaterData = permute(WaterData, [3 1 2]);

% MM (170505)
firstpoint_water = conj(WaterData(:,1,:));
channels_scale = squeeze(sqrt(sum(firstpoint_water .* conj(firstpoint_water),1)));
channels_scale = repmat(channels_scale, [1 nreceivers MRS_struct.p.npoints(ii)]);
channels_scale = permute(channels_scale, [2 3 1]);
firstpoint_water = repmat(firstpoint_water, [1 MRS_struct.p.npoints(ii) 1]) ./ channels_scale;

WaterData = WaterData .* firstpoint_water;
WaterData = squeeze(sum(WaterData,1));
MRS_struct.fids.data_water = WaterData;

% Use first point of water data to phase water-suppressed data
firstpoint = mean(firstpoint_water,3);
firstpoint = repmat(firstpoint, [1 1 size(FullData,3)]);

FullData = FullData .* firstpoint;
FullData = squeeze(sum(FullData,1));
MRS_struct.fids.data = FullData;

% Rescale, otherwise numbers blow up
MRS_struct.fids.data = MRS_struct.fids.data/1e11;
MRS_struct.fids.data_water = MRS_struct.fids.data_water/1e11;

end



