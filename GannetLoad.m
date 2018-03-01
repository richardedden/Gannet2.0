function MRS_struct = GannetLoad(metabfile, waterfile)
%Gannet 3.0 GannetLoad
%Started by RAEE Nov 5, 2012
%Updates by MGS, MM, GO 2016-2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Work flow summary
%   1. Pre-initialise
%   2. Determine data parameters from headers
%   3. Some housekeeping
%   4. Load data from files
%   5. Apply appropriate pre-processing
%   6. Output processed spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   0. Check the file list for typos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

missing=0;
for filecheck=1:length(metabfile)
    if ~exist(metabfile{filecheck},'file')
        disp(['The file ' metabfile{filecheck} ' (' num2str(filecheck) ')' ' is missing. Typo?'])
        missing=1;
    end
end
if nargin > 1
    for filecheck=1:length(waterfile)
        if ~exist(waterfile{filecheck},'file')
            disp(['The file ' waterfile(filecheck) ' is missing. Typo?'])
            missing=1;
        end
    end
end
if missing
    error('Not all the files are there, so I give up.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Pre-initialise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRS_struct.version.load = '171120'; % set to date when final updates have been made
MRS_struct.ii = 0;
MRS_struct.metabfile = metabfile;
MRS_struct = GannetPreInitialise(MRS_struct);

if MRS_struct.p.PRIAM % deciding how many voxels there are -- MGSaleh 2016
    vox = {MRS_struct.p.Vox};
else
    vox = {MRS_struct.p.Vox{1}};
end

if MRS_struct.p.HERMES % MGSaleh & MM 2016: for HERMES of GSH/Lac and GABAGlx/GSH
    % Swapping variables' values helps us with GannetLoad output -- MGSaleh 2016
    if strcmp(MRS_struct.p.target, 'Lac') && strcmp(MRS_struct.p.target2, 'GSH')
        [MRS_struct.p.target, MRS_struct.p.target2] = deal(MRS_struct.p.target2, MRS_struct.p.target);
    end
    if strcmp(MRS_struct.p.target, 'GSH') && strcmp(MRS_struct.p.target2, 'GABAGlx')
        [MRS_struct.p.target, MRS_struct.p.target2] = deal(MRS_struct.p.target2, MRS_struct.p.target);
    end
end

% Check whether or not there are water data
if nargin > 1
    MRS_struct.waterfile = waterfile;
    MRS_struct.p.Reference_compound = 'H2O';
else
    MRS_struct.p.Reference_compound = 'Cr';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2. Determine data parameters from header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iscell(metabfile) == 1 % it's a cell array, so work out the number of elements
    numpfiles = numel(metabfile);
    pfiles = metabfile;
else
    numpfiles = 1; % it's just one pfile
    pfiles{1} = metabfile;
end

MRS_struct = GannetDiscernDatatype(pfiles{1}, MRS_struct);

if strcmpi(MRS_struct.p.vendor,'Siemens')
    numpfiles = numpfiles/2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   3. Some housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create output folder
if ~exist('GannetLoad_output','dir')
    mkdir GannetLoad_output;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   4. Load data from files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:numpfiles % Loop over all files in the batch (from metabfile)
    
    MRS_struct.ii = ii;
    
    switch MRS_struct.p.vendor
        
        case 'GE'
            MRS_struct = GERead(MRS_struct, metabfile{ii});
            WaterData = MRS_struct.fids.data_water;
            MRS_struct.fids.data = MRS_struct.fids.data*MRS_struct.p.nrows(ii)/MRS_struct.p.Navg(ii);
            FullData = MRS_struct.fids.data;
            % Set up vector of which rows of data are ONs and OFFs
            switch MRS_struct.p.ONOFForder
                case 'onfirst'
                    if MRS_struct.p.HERMES % HERMES: GABAGlx or Lac and GSH (MM: 171120)
                        if strcmpi(MRS_struct.p.target, 'GABAGlx') && strcmpi(MRS_struct.p.target2, 'GSH')
                            % 1=ExpA, 2=ExpB, 3=ExpC, 4=ExpD (MM: 171120)
                            MRS_struct.fids.ON_OFF  = repmat([1 1 0 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                            MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                            MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                            MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                        end
                    else
                        MRS_struct.fids.ON_OFF = repmat([1 0], [1 size(MRS_struct.fids.data,2)/2]);
                    end
                case 'offfirst'
                    if MRS_struct.p.HERMES % HERMES: GABAGlx or Lac and GSH (MM: 171120)
                        if strcmpi(MRS_struct.p.target, 'GABAGlx') && strcmpi(MRS_struct.p.target2, 'GSH')
                            % 1=?, 2=?, 3=?, 4=? (MM: 171120)
                            MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                            MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                            MRS_struct.fids.ON_OFF  = repmat([1 0 0 1], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                            MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                        end
                    else
                        MRS_struct.fids.ON_OFF = repmat([0 1], [1 size(MRS_struct.fids.data,2)/2]);
                    end
            end
            
        case 'Siemens_twix'
            if exist('waterfile','var')
                MRS_struct = SiemensTwixRead(MRS_struct, metabfile{ii}, waterfile{ii});
                WaterData = MRS_struct.fids.data_water;
            else
                MRS_struct = SiemensTwixRead(MRS_struct, metabfile{ii});
            end
            % MM (160914): Need to set Water_Positive based on water signal
            if MRS_struct.p.Water_Positive == 0
                MRS_struct.fids.data = -MRS_struct.fids.data;
            end
            FullData = MRS_struct.fids.data;
            %Set up vector of which rows of data are ONs and OFFs
            switch MRS_struct.p.ONOFForder
                case 'onfirst'
                    MRS_struct.fids.ON_OFF = repmat([1 0],[1 size(MRS_struct.fids.data,2)/2]);
                case 'offfirst'
                    MRS_struct.fids.ON_OFF = repmat([0 1],[1 size(MRS_struct.fids.data,2)/2]);
            end
            
        case 'Siemens_dicom' % GO 11/01/2016
            if exist('waterfile','var')
                % Load the data. GannetLoad specifies a file in the first
                % place, so take apart that filename, and feed the containing
                % folder into the SiemensDICOMRead function. % GO 11/01/2016
                [imafolder,~,~] = fileparts(metabfile{ii}); % GO 02/05/2017
                [waterfolder,~,~] = fileparts(waterfile{ii}); % GO 02/05/2017
                MRS_struct = SiemensDICOMRead(MRS_struct,imafolder,waterfolder); % GO 02/05/2017
                MRS_struct.p.Reference_compound='H2O';
                WaterData = MRS_struct.fids.data_water;
            else
                MRS_struct.p.Reference_compound='Cr';
                % Same as above, but without parsing the waterfolder. % GO 02/05/2017
                [imafolder,~,~] = fileparts(metabfile{ii}); % GO 11/01/2016
                MRS_struct = SiemensDICOMRead(MRS_struct,imafolder); % GO 11/01/2016
            end
            FullData = MRS_struct.fids.data;
            
            % Fill up fields required for downstream processing % GO 11/01/2016
            switch MRS_struct.p.ONOFForder
                case 'onfirst'
                    MRS_struct.fids.ON_OFF=repmat([1 0],[1 MRS_struct.p.Navg(ii)/2]);
                    MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
                case 'offfirst'
                    MRS_struct.fids.ON_OFF=repmat([0 1],[1 MRS_struct.p.Navg(ii)/2]);
                    MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
            end
            
        case 'dicom' % GO 11/30/2016
            % care about water-unsuppressed files later % GO 11/30/2016
            if exist('waterfile','var')
                % Load the data. GannetLoad specifies a file in the first
                % place, so take apart that filename, and feed the containing
                % folder into the DICOMRead function. % GO 11/01/2016
                [dcmfolder,~,~] = fileparts(metabfile{ii}); % GO 02/05/2017
                [waterfolder,~,~] = fileparts(waterfile{ii}); % GO 02/05/2017
                MRS_struct = DICOMRead(MRS_struct,dcmfolder,waterfolder); % GO 02/05/2017
                MRS_struct.p.Reference_compound='H2O';
                WaterData = MRS_struct.fids.data_water;
            else
                MRS_struct.p.Reference_compound='Cr';
                % Same as above, but without parsing the waterfolder. % GO 02/05/2017
                [dcmfolder,~,~] = fileparts(metabfile{ii}); % GO 11/01/2016
                MRS_struct = DICOMRead(MRS_struct,dcmfolder); % GO 11/01/2016
            end
            FullData = MRS_struct.fids.data;
            
            % fill up fields required for downstream processing % GO 11/30/2016
            switch MRS_struct.p.ONOFForder
                % Not sure whether this is always the case, but the CMRR
                % sequence appears to go OFF-OFF-ON-ON in the DICOM
                % sorting?! Fixing this hard for now. GO 112017
                case 'onfirst'
                    if strcmp(MRS_struct.p.seq,'""%CustomerSeq%\eja_svs_mpress""')
                        MRS_struct.fids.ON_OFF=repmat([1 1 0 0],[1 MRS_struct.p.Navg(ii)/4]);
                        MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
                    else
                        MRS_struct.fids.ON_OFF=repmat([1 0],[1 MRS_struct.p.Navg(ii)/2]);
                        MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
                    end
                case 'offfirst'
                    if strcmp(MRS_struct.p.seq,'""%CustomerSeq%\eja_svs_mpress""')
                        MRS_struct.fids.ON_OFF=repmat([0 0 1 1],[1 MRS_struct.p.Navg(ii)/4]);
                        MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
                    else
                        MRS_struct.fids.ON_OFF=repmat([0 1],[1 MRS_struct.p.Navg(ii)/2]);
                        MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
                    end
            end
            
        case 'Siemens'
            if exist('waterfile','var')
                MRS_struct.p.Reference_compound = 'H2O';
                switch MRS_struct.p.ONOFForder
                    case 'offfirst'
                        MRS_struct = SiemensRead(MRS_struct, metabfile{ii*2-1},metabfile{ii*2}, waterfile{ii});
                    case 'onfirst'
                        MRS_struct = SiemensRead(MRS_struct, metabfile{ii*2},metabfile{ii*2-1}, waterfile{ii});
                end
                MRS_struct.p.Nwateravg = 1;
            else
                MRS_struct.p.Reference_compound = 'Cr';
                switch MRS_struct.p.ONOFForder
                    case 'offfirst'
                        MRS_struct = SiemensRead(MRS_struct, metabfile{ii*2-1},metabfile{ii*2});
                    case 'onfirst'
                        MRS_struct = SiemensRead(MRS_struct, metabfile{ii*2},metabfile{ii*2-1});
                end
            end
            FullData = MRS_struct.fids.data;
            if strcmp(MRS_struct.p.Reference_compound,'H2O')
                WaterData = MRS_struct.fids.data_water;
            end
            % Data are always read in OFF then ON
            switch MRS_struct.p.ONOFForder
                case 'onfirst'
                    MRS_struct.fids.ON_OFF = [1 0];
                case 'offfirst'
                    MRS_struct.fids.ON_OFF = [0 1];
            end
            
        case 'Philips'
            if strcmpi(MRS_struct.p.Reference_compound,'H2O')
                MRS_struct = PhilipsRead(MRS_struct, metabfile{ii}, waterfile{ii});
                WaterData = MRS_struct.fids.data_water;
            else
                MRS_struct = PhilipsRead(MRS_struct, metabfile{ii});
            end
            % Need to set Water_Positive based on water signal
            if MRS_struct.p.Water_Positive == 0
                MRS_struct.fids.data = -MRS_struct.fids.data;
            end
            FullData = MRS_struct.fids.data;
            if strcmpi(MRS_struct.p.seqorig,'JHU')
                switch MRS_struct.p.ONOFForder
                    case 'onfirst'
                        if MRS_struct.p.HERMES % HERMES: GABAGlx or Lac and GSH -- Added by MGSaleh & MM 2016
                            if strcmpi(MRS_struct.p.target, 'GABAGlx') && strcmpi(MRS_struct.p.target2, 'GSH')
                                % 1=?, 2=?, 3=?, 4=? (MM: 170703)
                                MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                                MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                            elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                                MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                                MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                            end
                        else
                            MRS_struct.fids.ON_OFF = repmat([1 0], [1 size(MRS_struct.fids.data,2)/2]);
                        end
                    case 'offfirst'
                        if MRS_struct.p.HERMES % HERMES: GABAGlx or Lac and GSH -- Added by MGSaleh & MM 2016
                            if strcmpi(MRS_struct.p.target, 'GABAGlx') && strcmpi(MRS_struct.p.target2, 'GSH')
                                % 1=ExpC, 2=ExpB, 3=ExpA, 4=ExpD (MM: 170703)
                                MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                                MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                            elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                                MRS_struct.fids.ON_OFF  = repmat([1 0 0 1], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                                MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                            end
                        else
                            MRS_struct.fids.ON_OFF = repmat([0 1], [1 size(MRS_struct.fids.data,2)/2]);
                        end
                end
            elseif strcmpi(MRS_struct.p.seqorig,'Philips')
                switch MRS_struct.p.ONOFForder
                    case 'onfirst'
                        MRS_struct.fids.ON_OFF = [ones(1,size(MRS_struct.fids.data,2)/2) zeros(1,size(MRS_struct.fids.data,2)/2)];
                    case 'offfirst'
                        MRS_struct.fids.ON_OFF = [zeros(1,size(MRS_struct.fids.data,2)/2) ones(1,size(MRS_struct.fids.data,2)/2)];
                end
            else 
                error('Philips sequence type not correctly specified in GannetPreInitialise. Options are JHU and Philips.') 
            end
            
        case 'Philips_data'
            if exist('waterfile','var')
                MRS_struct.p.Reference_compound = 'H2O';
                MRS_struct = PhilipsRead_data(MRS_struct, metabfile{ii}, waterfile{ii});
            else
                MRS_struct.p.Reference_compound = 'Cr';
                MRS_struct = PhilipsRead_data(MRS_struct, metabfile{ii});
            end
            if strcmpi(MRS_struct.p.Reference_compound,'H2O')
                WaterData = MRS_struct.fids.data_water;
            end
            FullData = MRS_struct.fids.data;
            switch MRS_struct.p.ONOFForder
                case 'onfirst'
                    MRS_struct.fids.ON_OFF = repmat([1 0],[MRS_struct.p.Navg(ii)/MRS_struct.p.nrows(ii) MRS_struct.p.nrows(ii)/2]);
                    MRS_struct.fids.ON_OFF = MRS_struct.fids.ON_OFF(:).';
                case 'offfirst'
                    MRS_struct.fids.ON_OFF = repmat([0 1],[MRS_struct.p.Navg(ii)/MRS_struct.p.nrows(ii) MRS_struct.p.nrows(ii)/2]);
                    MRS_struct.fids.ON_OFF = MRS_struct.fids.ON_OFF(:).';
            end
            
    end % end of vendor switch loop for data load
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   5. Apply appropriate pre-processing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % MM (160919): Zero-fill to obtain nominal spectral resolution of 0.061 Hz/point
    MRS_struct.p.ZeroFillTo(ii) = round(32768/2000*MRS_struct.p.sw(ii)); % MM (170727): round in case of non-integers
    MRS_struct.p.zf = MRS_struct.p.ZeroFillTo(ii)/MRS_struct.p.npoints(ii);
    time = (1:1:size(FullData,1))/MRS_struct.p.sw(ii);
    
    % Finish processing water data
    if strcmpi(MRS_struct.p.Reference_compound,'H2O')
        for kk = 1:length(vox)
            if strcmpi(MRS_struct.p.vendor,'GE')
                ComWater = mean(WaterData,2);
            elseif strcmpi(MRS_struct.p.vendor,'Siemens')
                ComWater = WaterData;
            elseif strcmpi(MRS_struct.p.vendor,'Siemens_twix')
                ComWater = WaterData;
            elseif (strcmpi(MRS_struct.p.vendor,'Siemens_dicom')) % GO 02/05/2017
                ComWater = mean(WaterData,2);
            elseif (strcmpi(MRS_struct.p.vendor,'dicom')) % GO 02/05/2017
                ComWater = mean(WaterData,2);
            elseif (strcmpi(MRS_struct.p.vendor,'Philips_raw')) % GO 02/05/2017
                ComWater = mean(WaterData(kk,:,:),2);
            else
                ComWater = WaterData.';
            end
            
            % Performing phase corrrection on the water-suppressed data
            % based on Klose (1990), MRM,14:26-30. The equation was
            % taken from Jiru (2008), EJR,67:202-217 -- MGSaleh 2016
            if MRS_struct.p.data_phase_correction
                if any(strcmpi(MRS_struct.p.vendor,{'Philips','Philips_data'}))
                    MRS_struct.fids.data = phase_correction_fids(MRS_struct.fids.data.', ComWater.');
                    MRS_struct.fids.data = MRS_struct.fids.data.';
                    FullData = MRS_struct.fids.data;
                else
                    MRS_struct.fids.data = phase_correction_fids(MRS_struct.fids.data.', ComWater);
                    MRS_struct.fids.data = MRS_struct.fids.data.';
                    FullData = MRS_struct.fids.data;
                end
            end
            
            % Performing phase corrrection on the unsuppressed water data
            if MRS_struct.p.water_phase_correction
                if any(strcmpi(MRS_struct.p.vendor,{'Philips','Philips_data'}))
                    ComWater = phase_correction_fids(ComWater.', ComWater.');
                    ComWater = ComWater.';
                else
                    ComWater = phase_correction_fids(ComWater, ComWater);
                end
            end
            
            % Line-broadening, zero-filling and FFT
            ComWater = ComWater .* exp(-time'*MRS_struct.p.LB*pi);
            MRS_struct.spec.(vox{kk}).water(ii,:) = fftshift(fft(ComWater,MRS_struct.p.ZeroFillTo(ii),1))';
        end
    end % end of H2O reference loop
    
    % Line-broadening, zero-filling and FFT
    FullData = FullData .* repmat((exp(-time'*MRS_struct.p.LB*pi)), [1 size(FullData,2)]);
    MRS_struct.fids.FullData = FullData;
    AllFramesFT = fftshift(fft(FullData,MRS_struct.p.ZeroFillTo(ii),1),1);
    
    % Work out frequency scale
    freqrange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
    MRS_struct.spec.freq = (MRS_struct.p.ZeroFillTo(ii)+1-(1:1:MRS_struct.p.ZeroFillTo(ii)))/MRS_struct.p.ZeroFillTo(ii)*freqrange+4.68-freqrange/2.0;
    % MM (170119)
    MRS_struct.p.df(ii) = abs(MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2));
    MRS_struct.p.SpecRes(ii) = MRS_struct.p.sw(ii)/MRS_struct.p.npoints(ii);
    MRS_struct.p.SpecResNominal(ii) = MRS_struct.p.sw(ii)/MRS_struct.p.ZeroFillTo(ii);
    MRS_struct.p.Tacq(ii) = 1/MRS_struct.p.SpecRes(ii);
    
    % Frame-by-frame determination of frequency of residual water (MM: 170201)
    water_range = MRS_struct.spec.freq-4.68 >= -0.2 & MRS_struct.spec.freq-4.68 <= 0.2;
    [~,FrameMaxPos] = max(abs(real(AllFramesFT(water_range,:))),[],1);
    % Not always true that water starts at 4.68, if drift is rapid...
    water_off = abs(MRS_struct.spec.freq-4.68);
    water_index = find(min(water_off)==water_off);
    % Determine frame shifts
    FrameShift = FrameMaxPos - water_index;
    
    switch MRS_struct.p.vendor
        case 'GE'
            AllFramesFTrealign = AllFramesFT;
        case {'Philips','Philips_data'}
            if any(strcmp(MRS_struct.p.target,{'GSH','GABAGlx'})) % Added by MGSaleh 2016
                AllFramesFTrealign = AllFramesFT;
            else
                for jj = 1:size(AllFramesFT,2)
                    AllFramesFTrealign(:,jj) = circshift(AllFramesFT(:,jj), -FrameShift(jj)); % is this used????
                end
            end
            %This quite possibly doesn't carry through, as it seems
            %that the later stuff all starts with AllFramesFT, not
            %AllFramesFTrealign
        case 'Siemens_twix'
            AllFramesFTrealign = AllFramesFT;
        case 'Siemens_dicom'
            AllFramesFTrealign = AllFramesFT;
        case 'dicom'
            AllFramesFTrealign = AllFramesFT;
    end
    
    % MM (170703)
    freqWaterRange = MRS_struct.spec.freq(water_range);
    MRS_struct.fids.waterfreq(ii,:) = freqWaterRange(FrameMaxPos);
    
    % MM (170629): Estimate average amount of F0 offset
    if any(strcmp(MRS_struct.p.vendor,{'Siemens','Siemens_twix','Siemens_dicom'}))
        MRS_struct.out.AvgDeltaF0(ii) = mean(freqWaterRange(FrameMaxPos) - 4.7); % Siemens assumes 4.7 ppm as F0
    else
        MRS_struct.out.AvgDeltaF0(ii) = mean(freqWaterRange(FrameMaxPos) - 4.68);
    end
    
    % Frame-by-frame alignment
    switch MRS_struct.p.AlignTo
        case 'Cr'
            [AllFramesFTrealign, MRS_struct] = AlignUsingPeak(AllFramesFTrealign,MRS_struct);
            %AllFramesFTrealign = AlignUsingCr(AllFramesFTrealign,MRS_struct.p.ONOFForder,n);
        case 'Cho'
            [AllFramesFTrealign, MRS_struct] = AlignUsingPeak(AllFramesFTrealign,MRS_struct);
        case 'H2O'
            [AllFramesFTrealign, MRS_struct] = AlignUsingH2O(AllFramesFTrealign,MRS_struct);
        case 'NAA'
            [AllFramesFTrealign, MRS_struct] = AlignUsingPeak(AllFramesFTrealign,MRS_struct);
        case 'SpecReg'
            [AllFramesFTrealign, MRS_struct] = Spectral_Registration(MRS_struct,0);
        case 'SpecRegDual'
            %Dual-channel Spectral Registration is applied separately to ON and OFF and they are coregistered after...
            [AllFramesFTrealign, MRS_struct] = Spectral_Registration(MRS_struct,0,1);
        case 'SpecRegHERMES' % MM (170703)
            [AllFramesFTrealign, MRS_struct] = Spectral_Registration_HERMES(MRS_struct);
    end % end of switch for alignment
    
    % Separate ON/OFF data and generate DIFF spectra
    for kk = 1:length(vox) % loop over voxels -- MGSaleh 2016
        if MRS_struct.p.HERMES % MGSaleh 2016, MM (170703)
            
            % Target 1: GABA or GSH
            OFF = mean(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0)' & MRS_struct.out.reject(:,ii)==0), 2);
            ON  = mean(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1)' & MRS_struct.out.reject(:,ii)==0), 2);
            
            MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).off(ii,:) = OFF;
            MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).on(ii,:)  = ON;
            
            MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) = (ON-OFF)/2;
            MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:) = (mean(AllFramesFT(:,MRS_struct.fids.ON_OFF==1),2) - mean(AllFramesFT(:,MRS_struct.fids.ON_OFF==0),2))/2;
            
            % Target 2: GSH or Lac
            OFF2 = mean(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF2==0)' & MRS_struct.out.reject(:,ii)==0), 2);
            ON2  = mean(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF2==1)' & MRS_struct.out.reject(:,ii)==0), 2);
            
            MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).off(ii,:) = OFF2;
            MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).on(ii,:)  = ON2;
            
            MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff(ii,:) = (ON2-OFF2)/2;
            MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff_noalign(ii,:) = (mean(AllFramesFT(:,MRS_struct.fids.ON_OFF2==1),2) - mean(AllFramesFT(:,MRS_struct.fids.ON_OFF2==0),2))/2;
            
            % Remove residual water from diff and diff_noalign spectra using HSVD -- GO & MGSaleh 2016
            if MRS_struct.p.water_removal

                % Convert DIFF spectra to time domain, apply water filter, convert back to frequency domain
                MRS_struct.fids.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) = waterremovalSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:).')), ...
                    MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048); % MM (171121)
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) = fftshift(fft(MRS_struct.fids.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:)));
                
                MRS_struct.fids.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff(ii,:) = waterremovalSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff(ii,:).')), ...
                    MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048);
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff(ii,:) = fftshift(fft(MRS_struct.fids.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff(ii,:)));
                
                MRS_struct.fids.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:) = waterremovalSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:).')), ...
                    MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048);
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:) = fftshift(fft(MRS_struct.fids.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:)));
                
                MRS_struct.fids.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff_noalign(ii,:) = waterremovalSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff_noalign(ii,:).')), ...
                    MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048);
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff_noalign(ii,:) = fftshift(fft(MRS_struct.fids.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff_noalign(ii,:)));
                
                % MM (170703): Need to perform baseline correction on filtered data
                freqbounds = MRS_struct.spec.freq <= 10 & MRS_struct.spec.freq >= 9;
                baseMean_diff1 = mean(real(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,freqbounds)));
                baseMean_diffnoalign1 = mean(real(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,freqbounds)));
                baseMean_diff2 = mean(real(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff(ii,freqbounds)));
                baseMean_diffnoalign2 = mean(real(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff_noalign(ii,freqbounds)));
                
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) = MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) - baseMean_diff1;
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:) = MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:) - baseMean_diffnoalign1;
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff(ii,:) = MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff(ii,:) - baseMean_diff2;
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff_noalign(ii,:) = MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).diff_noalign(ii,:) - baseMean_diffnoalign2;
            end
            
        else
            
            if strcmp(MRS_struct.p.target, 'GSH')
                
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).off(ii,:) = mean(AllFramesFTrealign(:,((MRS_struct.fids.ON_OFF==0)'&(MRS_struct.out.reject(:,ii)==0))),2);
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).on(ii,:) = mean(AllFramesFTrealign(:,((MRS_struct.fids.ON_OFF==1)'&(MRS_struct.out.reject(:,ii)==0))),2);
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) = (MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).on(ii,:)-MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).off(ii,:))/2;
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:) = (mean(AllFramesFT(:,(MRS_struct.fids.ON_OFF==1)),2)-mean(AllFramesFT(:,(MRS_struct.fids.ON_OFF==0)),2))/2;
                
                %For GSH data, the residual water signal in the DIFF spectrum is
                %helpful for an additional phasing step... and messes up fitting
                %otherwise. MGSaleh 2016 moved it to this place for
                %completeness
                if ~strcmp(MRS_struct.p.AlignTo, 'SpecRegHERMES') % Don't apply residual phasing if SpecRegHERMES is used (MM: 180108)
                    residual_phase = pi-atan2(imag(sum(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:))),real(sum(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:))));
                    MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) = (MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:))*exp(1i*residual_phase);
                end
                
                if MRS_struct.p.Water_Positive == 0
                    MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) = -MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:);
                end
                
                if MRS_struct.p.water_removal
                    % Convert DIFF spectra to time domain, apply water filter, convert back to frequency domain
                    MRS_struct.fids.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) = waterremovalSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:).')), ...
                        MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048); % MM (171121)
                    MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) = fftshift(fft(MRS_struct.fids.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:)));
                    
                    MRS_struct.fids.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:) = waterremovalSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:).')), ...
                        MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048);
                    MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:) = fftshift(fft(MRS_struct.fids.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:)));
                    
                    % MM (170703): Need to perform baseline correction on filtered data
                    freqbounds = MRS_struct.spec.freq <= 10 & MRS_struct.spec.freq >= 9;
                    baseMean_diff1 = mean(real(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,freqbounds)));
                    baseMean_diffnoalign1 = mean(real(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,freqbounds)));
                    
                    MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) = MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) - baseMean_diff1;
                    MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:) = MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:) - baseMean_diffnoalign1;
                    
                end
            else
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).off(ii,:) = mean(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0)' & MRS_struct.out.reject(:,ii)==0), 2);
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).on(ii,:)  = mean(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1)' & MRS_struct.out.reject(:,ii)==0), 2);
                
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:) = (MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).on(ii,:) - MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).off(ii,:))/2;
                MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff_noalign(ii,:) = (mean(AllFramesFT(:,(MRS_struct.fids.ON_OFF==1)),2) - mean(AllFramesFT(:,(MRS_struct.fids.ON_OFF==0)),2))/2;
                
            end
            
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   6. Build GannetLoad Output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ishandle(101)
        clf(101); % MM (170629)
    end
    h = figure(101);
    % MM (170629): Open figure in center of screen
    scr_sz = get(0, 'ScreenSize');
    fig_w = 1000;
    fig_h = 707;
    set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
    set(h,'Color',[1 1 1]);
    figTitle = 'GannetLoad Output';
    set(gcf,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
    
    % Top left
    ha = subplot(2,2,1);
    GannetPlotPrePostAlign(MRS_struct, vox, ii);
    if MRS_struct.p.HERMES
        title({'Edited Spectra';'(pre- and post-alignment)'});
    else
        title({'Edited Spectrum';'(pre- and post-alignment)'});
    end
    xlabel('ppm');
    set(gca,'YTick',[]);
    
    % Top right
    hb = subplot(2,2,2);
    rejectframesplot = (1./MRS_struct.out.reject(:,ii).') .*  MRS_struct.fids.waterfreq(ii,:);
    plot(1:size(FullData,2), MRS_struct.fids.waterfreq(ii,:)', '-', 1:size(FullData,2), rejectframesplot, 'ro');
    set(gca,'XLim',[0 size(FullData,2)]);
    xlabel('average'); ylabel('\omega_0');
    title('Water Frequency, ppm');
    
    % Bottom left
    hc = subplot(2,2,3);
    if ~strcmp(MRS_struct.p.AlignTo,'no')
        CrFitLimLow = 2.72;
        CrFitLimHigh = 3.12;
        plotrange = MRS_struct.spec.freq <= CrFitLimHigh & MRS_struct.spec.freq >= CrFitLimLow; % MM (170705)
        CrFitRange = sum(plotrange);
        plotrealign = [real(AllFramesFT(plotrange,:)); real(AllFramesFTrealign(plotrange,:))];
        % Don't display rejects
        plotrealign(CrFitRange+1:end,(MRS_struct.out.reject(:,ii).'==1))=min(plotrealign(:));
        imagesc(plotrealign);
        title({'Cr Frequency','(pre- and post-alignment)'});
        xlabel('average');
        set(gca,'YTick', [1 CrFitRange CrFitRange+CrFitRange*(CrFitLimHigh-3.02)/(CrFitLimHigh-CrFitLimLow) CrFitRange*2]);
        set(gca,'YTickLabel', [CrFitLimHigh CrFitLimLow 3.02 CrFitLimLow]);
        % Add in labels for pre/post
        text(size(plotrealign,2)/18*17,0.4*size(plotrealign,1), 'PRE', 'Color', [1 1 1], 'HorizontalAlignment', 'right');
        text(size(plotrealign,2)/18*17,0.9*size(plotrealign,1), 'POST', 'Color', [1 1 1], 'HorizontalAlignment', 'right');
    else
        tmp = 'No realignment';
        text(0, 0.9, tmp, 'FontName', 'Courier');
    end
    
    % Bottom right
    hd = subplot(2,2,4);
    axis off;
    
    text_pos = 0.9;
    
    % MM (180112)
    if strcmp(MRS_struct.p.vendor,'Siemens')
        [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii*2-1});
    else
        [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii});
    end
    
    text(0, text_pos, 'Filename', 'FontName', 'Helvetica', 'FontSize', 13);
    text(0.275, text_pos, [': ' tmp tmp2], 'FontName', 'Helvetica', 'FontSize', 13, 'Interpreter', 'none');
    
    tmp = [': ' num2str(MRS_struct.p.Navg(ii)) ' averages'];
    text(0, text_pos - 0.1, 'Navg', 'FontName', 'Helvetica', 'FontSize', 13);
    text(0.275, text_pos - 0.1, tmp, 'FontName', 'Helvetica', 'FontSize', 13);
    
    if isfield(MRS_struct.p,'voxdim')
        tmp = [': '  num2str(MRS_struct.p.voxdim(ii,1)*MRS_struct.p.voxdim(ii,2)*MRS_struct.p.voxdim(ii,3)/1e3) ' mL'];
        text(0, text_pos - 0.2, 'Volume', 'FontName', 'Helvetica', 'FontSize', 13);
        text(0.275, text_pos - 0.2, tmp, 'FontName', 'Helvetica', 'FontSize', 13);
    end
    
    tmp = [': '  MRS_struct.p.AlignTo];
    text(0, text_pos - 0.3, 'Alignment', 'FontName', 'Helvetica', 'FontSize', 13);
    text(0.275, text_pos - 0.3, tmp, 'FontName', 'Helvetica', 'FontSize', 13);
    
    tmp = [': ' num2str(MRS_struct.p.LB,2) ' Hz'];
    text(0, text_pos - 0.4, 'LB', 'FontName', 'Helvetica', 'FontSize', 13);
    text(0.275, text_pos - 0.4, tmp, 'FontName', 'Helvetica', 'FontSize', 13);
    
    tmp = [': '  num2str(sum(MRS_struct.out.reject(:,ii),1)) ];
    text(0, text_pos - 0.5, 'Rejects', 'FontName', 'Helvetica', 'FontSize', 13);
    text(0.275, text_pos - 0.5, tmp, 'FontName', 'Helvetica', 'FontSize', 13);
    
    tmp = [': ' MRS_struct.version.load];
    text(0, text_pos - 0.6, 'LoadVer', 'FontName', 'Helvetica', 'FontSize', 13);
    text(0.275, text_pos - 0.6, tmp, 'FontName', 'Helvetica', 'FontSize', 13);
    
    % Add Gannet logo
    Gannet_path = which('GannetLoad');
    Gannet_logo = [Gannet_path(1:end-13) '/Gannet3_logo.png'];
    A2 = imread(Gannet_logo,'png','BackgroundColor',[1 1 1]);
    axes('Position',[0.80, 0.05, 0.15, 0.15]);
    image(A2);
    axis off;
    axis square;
        
    % For Philips .data
    if strcmpi(MRS_struct.p.vendor,'Philips_data')
        fullpath = MRS_struct.metabfile{ii};
        fullpath = regexprep(fullpath, '.data', '_data');
        fullpath = regexprep(fullpath, '\', '_');
        fullpath = regexprep(fullpath, '/', '_');
    end
        
    % MM (180112)
    if strcmp(MRS_struct.p.vendor,'Siemens')
        [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii*2-1});
    else
        [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii});
    end
    
    if sum(strcmp(listfonts,'Helvetica')) > 0
        set([ha,hb,hc,hd],'FontName','Helvetica'); % GO 11/16/2017; MM: 171120
    end
    
    % Save PDF
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize',[11 8.5]);
    set(gcf,'PaperPosition',[0 0 11 8.5]);
    
    if strcmpi(MRS_struct.p.vendor,'Philips_data')
        pdfname = fullfile('GannetLoad_output', [fullpath '_load.pdf']); % MM (180112)
    else
        pdfname = fullfile('GannetLoad_output', [metabfile_nopath '_load.pdf']); % MM (180112)
    end
    saveas(h, pdfname);
    
    % Save the processed data into an SDAT file
    if MRS_struct.p.sdat
        if strcmpi(MRS_struct.p.vendor,'Philips')
            if strcmpi(MRS_struct.p.vendor,'Philips_data')
                %sdat_G_name=[ 'MRSload_output/' fullpath '_G.data' ]
                %NOT SUPPORTED
            else
                %set up filenames for sdat output
                sdat_G_name=['GannetLoad_output/' metabfile_nopath  '_G.sdat'];
                spar_G_name=['GannetLoad_output/' metabfile_nopath  '_G.spar'];
                %make file copies for sdat output
                copyfile(metabfile{ii},sdat_G_name);
                sparname=metabfile{ii};
                sparname = [sparname(1:(end-4)) MRS_struct.p.spar_string];
                copyfile(sparname,spar_G_name);
                %write into the sdat file
                %What do we write
                sdat_diff_out=conj(ifft(fftshift(MRS_struct.spec.diff(ii,:),2),[],2));
                sdat_diff_out=sdat_diff_out(1:MRS_struct.p.npoints(ii));
                %Also write out OFF
                sdat_off_out=conj(ifft(fftshift(MRS_struct.spec.GABA.off(ii,:),2),[],2));
                sdat_off_out=sdat_off_out(1:MRS_struct.p.npoints(ii));
                %How do we write it out?
                fileid  = fopen(sdat_G_name,'w','ieee-le');
                ff(:,1:2:2*MRS_struct.p.npoints(ii)) = real(sdat_diff_out);
                ff(:,2:2:2*MRS_struct.p.npoints(ii)) = imag(sdat_diff_out);
                gg(:,1:2:2*MRS_struct.p.npoints(ii)) = real(sdat_off_out);
                gg(:,2:2:2*MRS_struct.p.npoints(ii)) = imag(sdat_off_out);
                fwriteVAXD(fileid,[ff.' gg.'],'float');
                fclose(fileid);
            end
        end
    end
    
    % 140116: ADH reorder structure
    if(isfield(MRS_struct, 'waterfile') == 1)
        structorder = {'version', 'ii', 'metabfile', ...
            'waterfile', 'p', 'fids', 'spec', 'out'};
    else
        structorder = {'version', 'ii', 'metabfile', ...
            'p', 'fids', 'spec', 'out'};
    end
    MRS_struct = orderfields(MRS_struct, structorder);
    
end % end of load-and-processing loop over datasets

end


