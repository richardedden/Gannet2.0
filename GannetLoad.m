function MRS_struct = GannetLoad(metabfile, waterfile)
% Gannet 3.0 GannetLoad
% Started by RAEE Nov 5, 2012
% Updates by MGS, MM, GO 2016-2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Workflow summary
%   1. Pre-initialise
%   2. Determine data parameters from headers
%   3. Load data from files
%   4. Reconstruction of coil-sensitivity maps (PRIAM only)
%   5. Apply appropriate pre-processing
%   6. Output processed spectra
%   7. Build GannetLoad output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   0. Check the file list for typos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

missing = 0;
for filecheck = 1:length(metabfile)
    % If only water-suppressed data are provided, select Cr as reference.
    MRS_struct.p.Reference_compound = 'Cr';
    if ~exist(metabfile{filecheck},'file')
        disp(['The file ' metabfile{filecheck} ' (' num2str(filecheck) ')' ' is missing. Typo?']);
        missing = 1;
    end
end
if nargin > 1
    % If water-unsuppressed data are provided, select H2O as reference.
    MRS_struct.waterfile = waterfile;
    MRS_struct.p.Reference_compound = 'H2O';
    for filecheck=1:length(waterfile)
        if ~exist(waterfile{filecheck},'file')
            disp(['The file ' waterfile(filecheck) ' is missing. Typo?'])
            missing = 1;
        end
    end
end
if missing
    error('Not all the files are there, so I give up.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Pre-initialise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRS_struct.version.load = '180326'; % set to date when final updates have been made
MRS_struct.ii = 0;
MRS_struct.metabfile = metabfile;
MRS_struct = GannetPreInitialise(MRS_struct);

if MRS_struct.p.PRIAM % deciding how many voxels there are -- MGSaleh 2016
    vox = MRS_struct.p.Vox;
else
    vox = {MRS_struct.p.Vox{1}};
end

if MRS_struct.p.HERMES % MGSaleh & MM 2016: for HERMES of GSH/Lac and GABAGlx/GSH
    % Swapping variables' values helps us with GannetLoad output -- MGSaleh 2016
    if strcmp(MRS_struct.p.target, 'Lac') && strcmp(MRS_struct.p.target2, 'GSH')
        [MRS_struct.p.target, MRS_struct.p.target2] = deal(MRS_struct.p.target2, MRS_struct.p.target);
    end
    if strcmp(MRS_struct.p.target, 'GSH') && any(strcmp(MRS_struct.p.target2, {'GABA','Glx','GABAGlx'}))
        [MRS_struct.p.target, MRS_struct.p.target2] = deal(MRS_struct.p.target2, MRS_struct.p.target);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2. Determine data parameters from header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine number of provided water-suppressed files in the batch
MRS_struct.p.Reference_compound='Cr';
numscans = numel(metabfile);

% Discern input data format
MRS_struct = GannetDiscernDatatype(metabfile{1}, MRS_struct);

% For Siemens RDA, each acquisition has two RDA files, i.e. correct the
% number:
if strcmpi(MRS_struct.p.vendor,'Siemens_rda')
    numscans = numscans/2;
end
% Determine number of provided water-unsuppressed files in the batch
if exist('waterfile','var')
    MRS_struct.p.Reference_compound='H2O';
    numwaterscans = numel(waterfile);
    if numwaterscans ~= numscans
        error ('Number of water-unsuppressed files does not match number of water-suppressed files.');
    end
end

% Create output folder
if ~exist('GannetLoad_output','dir')
    mkdir GannetLoad_output;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   3. Load data from files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:numscans % Loop over all files in the batch (from metabfile)
    
    MRS_struct.ii = ii;
    
    switch MRS_struct.p.vendor
        
        case 'GE'
            MRS_struct = GERead(MRS_struct, metabfile{ii});
            WaterData = MRS_struct.fids.data_water;
            MRS_struct.p.Reference_compound='H2O';
            MRS_struct.fids.data = MRS_struct.fids.data*MRS_struct.p.nrows(ii)/MRS_struct.p.Navg(ii);
            FullData = MRS_struct.fids.data;
            % Determine order of ON and OFF acquisitions
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
            
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
            % Determine order of ON and OFF acquisitions
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
            FullData = MRS_struct.fids.data;
            
        case 'Siemens_dicom' % GO 11/01/2016
            if exist('waterfile','var')
                % Load the data. GannetLoad specifies a file in the first
                % place, so take apart that filename, and feed the containing
                % folder into the SiemensDICOMRead function. % GO 11/01/2016
                [imafolder,~,~] = fileparts(metabfile{ii}); % GO 02/05/2017
                [waterfolder,~,~] = fileparts(waterfile{ii}); % GO 02/05/2017
                MRS_struct = SiemensDICOMRead(MRS_struct,imafolder,waterfolder); % GO 02/05/2017
                WaterData = MRS_struct.fids.data_water;
            else
                % Same as above, but without parsing the waterfolder. % GO 02/05/2017
                [imafolder,~,~] = fileparts(metabfile{ii}); % GO 11/01/2016
                MRS_struct = SiemensDICOMRead(MRS_struct,imafolder); % GO 11/01/2016
            end
            FullData = MRS_struct.fids.data;
            % Determine order of ON and OFF acquisitions
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
            
        case 'dicom' % GO 11/30/2016
            if exist('waterfile','var')
                % Load the data. GannetLoad specifies a file in the first
                % place, so take apart that filename, and feed the containing
                % folder into the DICOMRead function. % GO 11/01/2016
                [dcmfolder,~,~] = fileparts(metabfile{ii}); % GO 02/05/2017
                [waterfolder,~,~] = fileparts(waterfile{ii}); % GO 02/05/2017
                MRS_struct = DICOMRead(MRS_struct,dcmfolder,waterfolder); % GO 02/05/2017
                WaterData = MRS_struct.fids.data_water;
            else
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
%                     if strcmp(MRS_struct.p.seq,'""%CustomerSeq%\eja_svs_mpress""')
%                         MRS_struct.fids.ON_OFF=repmat([1 1 0 0],[1 MRS_struct.p.Navg(ii)/4]);
%                         MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
%                     else
                        MRS_struct.fids.ON_OFF=repmat([1 0],[1 MRS_struct.p.Navg(ii)/2]);
                        MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
%                    end
                case 'offfirst'
%                     if strcmp(MRS_struct.p.seq,'""%CustomerSeq%\eja_svs_mpress""')
%                         MRS_struct.fids.ON_OFF=repmat([0 0 1 1],[1 MRS_struct.p.Navg(ii)/4]);
%                         MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
%                     else
                        MRS_struct.fids.ON_OFF=repmat([0 1],[1 MRS_struct.p.Navg(ii)/2]);
                        MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
%                    end
            end
            
        case 'Siemens_rda'
            if exist('waterfile','var')
                MRS_struct = SiemensRead(MRS_struct, metabfile{ii*2}, metabfile{ii*2-1}, waterfile{ii});
                WaterData = MRS_struct.fids.data_water;
                MRS_struct.p.Nwateravg = 1;
            else
                MRS_struct = SiemensRead(MRS_struct, metabfile{ii*2}, metabfile{ii*2-1});
            end
            FullData = MRS_struct.fids.data;
            % Determine order of ON and OFF acquisitions
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
            
        case 'Philips'
            if exist('waterfile','var')
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
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
            
        case 'Philips_data'
            % If a water reference scan is acquired, it is saved as a mix
            % in the DATA/LIST files. Later: add option to provide an additional
            % water reference file (i.e. short-TE). GO 03/02/2018
            if nargin > 1
                MRS_struct = PhilipsRead_data(MRS_struct, metabfile{ii},waterfile{ii});
            else
                MRS_struct = PhilipsRead_data(MRS_struct, metabfile{ii});
            end
            if isfield(MRS_struct.fids, 'data_water')
                MRS_struct.p.Reference_compound = 'H2O';
                WaterData = MRS_struct.fids.data_water;
            else
                MRS_struct.p.Reference_compound = 'Cr';
            end
            MRS_struct = SpecifyOnOffOrder(MRS_struct); %For 3T and 7T -- 08212018 MGSaleh
            FullData = MRS_struct.fids.data;
        
        case 'Philips_raw' % GO 11/01/2016
            
            MRS_struct = PhilipsRawLoad(MRS_struct,metabfile{ii},3,0); % GO 11/02/2016 
            MRS_struct.fids.data=conj(squeeze(MRS_struct.multivoxel.allsignals(:,:,1,:)));
            if exist('waterfile','var')
                MRS_struct.p.Reference_compound = 'H2O';
                WaterData = MRS_struct.fids.data_water; % GO 11/03/2016
            end % GO 11/03/2016
            FullData = MRS_struct.fids.data;
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
    
    end % end of vendor switch loop for data load
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   4. Reconstruction of coil-sensitivity maps
    %      (PRIAM only)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if a PRIAM dataset is processed, load the coil reference scan and
    % calculate the SENSE reconstruction matrix here
    if MRS_struct.p.PRIAM
        MRS_struct = senseRecon(MRS_struct);
        PRIAMData = zeros(length(MRS_struct.p.Vox),MRS_struct.p.Navg,MRS_struct.p.npoints);
        PRIAMWaterData = zeros(length(MRS_struct.p.Vox),MRS_struct.p.Nwateravg,MRS_struct.p.npoints);
        for kk = 1:MRS_struct.p.Navg
            PRIAMData(:,kk,:) = MRS_struct.p.SENSE.U * squeeze(FullData(:,kk,:));
            % Phase by multiplying with normalized complex conjugate of first point
            conj_norm = conj(PRIAMData(:,kk,1)) ./ abs(conj(PRIAMData(:,kk,1)));
            PRIAMData(:,kk,:) = PRIAMData(:,kk,:) .* repmat(conj_norm, [1 1 MRS_struct.p.npoints]);
        end
        for kk = 1:MRS_struct.p.Nwateravg
            PRIAMWaterData(:,kk,:) = MRS_struct.p.SENSE.U * squeeze(WaterData(:,kk,:));
            % Phase by multiplying with normalized complex conjugate of first point
            conj_norm = conj(PRIAMWaterData(:,kk,1)) ./ abs(conj(PRIAMWaterData(:,kk,1)));
            PRIAMWaterData(:,kk,:) = PRIAMWaterData(:,kk,:) .* repmat(conj_norm, [1 1 MRS_struct.p.npoints]);
        end
    elseif strcmp(MRS_struct.p.vendor,'Philips_data') && ceil(MRS_struct.p.LarmorFreq(MRS_struct.ii)) < 290 %Added info about averaging in the coil dimensions -- 08232018 MGSaleh
        FullData = squeeze(sum(FullData,1))';
        MRS_struct.fids.data = FullData;
        if strcmp(MRS_struct.p.Reference_compound,'H2O')
            WaterData = squeeze(sum(WaterData,1))';
            MRS_struct.fids.data_water = WaterData;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   5. Apply appropriate pre-processing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for kk = 1:length(vox) % loop over number of voxels % GO 03/26/2018
        % Select data from first voxel
        if MRS_struct.p.PRIAM
            FullData = squeeze(-PRIAMData(kk,:,:))';
            MRS_struct.fids.data = FullData;
            WaterData = squeeze(PRIAMWaterData(kk,:,:))';
            MRS_struct.fids.data_water = WaterData;
        end
        % MM (160919): Zero-fill to obtain nominal spectral resolution of 0.061 Hz/point
        MRS_struct.p.ZeroFillTo(ii) = round(32768/2000*MRS_struct.p.sw(ii)); % MM (170727): round in case of non-integers
        MRS_struct.p.zf = MRS_struct.p.ZeroFillTo(ii)/MRS_struct.p.npoints(ii);
        time = (1:1:size(FullData,1))/MRS_struct.p.sw(ii);
        
        % Finish processing water data
        if strcmpi(MRS_struct.p.Reference_compound,'H2O')
            if strcmpi(MRS_struct.p.vendor,'GE')
                ComWater = mean(WaterData,2);
            elseif strcmpi(MRS_struct.p.vendor,'Siemens_rda')
                ComWater = WaterData;
            elseif strcmpi(MRS_struct.p.vendor,'Siemens_twix')
                ComWater = WaterData;
            elseif (strcmpi(MRS_struct.p.vendor,'Siemens_dicom')) % GO 02/05/2017
                ComWater = mean(WaterData,2);
            elseif (strcmpi(MRS_struct.p.vendor,'dicom')) % GO 02/05/2017
                ComWater = mean(WaterData,2);
            elseif (strcmpi(MRS_struct.p.vendor,'Philips_raw')) % GO 02/05/2017
                ComWater = mean(WaterData(kk,:,:),2);
            elseif (strcmpi(MRS_struct.p.vendor,'Philips_data')) % GO 03/18/2018
                ComWater = mean(WaterData,2);
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
            % GO (180514): Water data may have different bandwidth
            if isfield(MRS_struct.p,'sw_water')
                time_water = (1:1:size(ComWater,1))/MRS_struct.p.sw_water(ii);
            else
                time_water = (1:1:size(ComWater,1))/MRS_struct.p.sw(ii);
            end
            ComWater = ComWater .* exp(-time_water'*MRS_struct.p.LB*pi);
            MRS_struct.spec.(vox{kk}).water(ii,:) = fftshift(fft(ComWater,MRS_struct.p.ZeroFillTo(ii),1))';
        end % end of H2O reference loop
        
        % Line-broadening, zero-filling and FFT
        FullData = FullData .* repmat((exp(-time'*MRS_struct.p.LB*pi)), [1 size(FullData,2)]);
        MRS_struct.fids.FullData = FullData;
        AllFramesFT = fftshift(fft(FullData,MRS_struct.p.ZeroFillTo(ii),1),1);
        
        % Work out frequency scale
        freqrange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
        if MRS_struct.p.phantom
            F0 = 4.8;
        else
            F0 = 4.68;
        end
        MRS_struct.spec.freq = (MRS_struct.p.ZeroFillTo(ii)+1-(1:1:MRS_struct.p.ZeroFillTo(ii)))/MRS_struct.p.ZeroFillTo(ii)*freqrange+F0-freqrange/2;
        % MM (170119)
        MRS_struct.p.df(ii) = abs(MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2));
        MRS_struct.p.SpecRes(ii) = MRS_struct.p.sw(ii)/MRS_struct.p.npoints(ii);
        MRS_struct.p.SpecResNominal(ii) = MRS_struct.p.sw(ii)/MRS_struct.p.ZeroFillTo(ii);
        MRS_struct.p.Tacq(ii) = 1/MRS_struct.p.SpecRes(ii);
        
        % Frame-by-frame determination of frequency of residual water (if MEGA-PRESS) or Cr (if HERMES) (MM: 180801)
        if MRS_struct.p.HERMES
            F0freqRange = MRS_struct.spec.freq - 3.02 >= -0.1 & MRS_struct.spec.freq - 3.02 <= 0.1;
        else
            F0freqRange = MRS_struct.spec.freq - F0 >= -0.2 & MRS_struct.spec.freq - F0 <= 0.2;
        end
        [~,FrameMaxPos] = max(abs(real(AllFramesFT(F0freqRange,:))),[],1);
        F0freqRange = MRS_struct.spec.freq(F0freqRange);
        MRS_struct.spec.F0freq(ii,:) = F0freqRange(FrameMaxPos);
        
        % MM (180801): Estimate average amount of F0 offset
        if MRS_struct.p.HERMES
            MRS_struct.out.AvgDeltaF0(ii) = mean(F0freqRange(FrameMaxPos) - 3.02);
        elseif any(strcmp(MRS_struct.p.vendor,{'Siemens_rda','Siemens_twix','Siemens_dicom'}))
            MRS_struct.out.AvgDeltaF0(ii) = mean(F0freqRange(FrameMaxPos) - 4.7); % Siemens assumes 4.7 ppm as F0
        else
            MRS_struct.out.AvgDeltaF0(ii) = mean(F0freqRange(FrameMaxPos) - F0);
        end
        
        AllFramesFTrealign = AllFramesFT;
        
        % Frame-by-frame alignment
        switch MRS_struct.p.AlignTo
            case 'Cr'
                [AllFramesFTrealign, MRS_struct] = AlignUsingPeak(AllFramesFTrealign,MRS_struct);
                %AllFramesFTrealign = AlignUsingCr(AllFramesFTrealign,MRS_struct.p.ONOFForder,n);
                MRS_struct.spec.AllFramesFTrealign = AllFramesFTrealign;
            case 'Cho'
                [AllFramesFTrealign, MRS_struct] = AlignUsingPeak(AllFramesFTrealign,MRS_struct);
                MRS_struct.spec.AllFramesFTrealign = AllFramesFTrealign;
            case 'H2O'
                [AllFramesFTrealign, MRS_struct] = AlignUsingH2O(AllFramesFTrealign,MRS_struct);
                MRS_struct.spec.AllFramesFTrealign = AllFramesFTrealign;
            case 'NAA'
                [AllFramesFTrealign, MRS_struct] = AlignUsingPeak(AllFramesFTrealign,MRS_struct);
                MRS_struct.spec.AllFramesFTrealign = AllFramesFTrealign;
            case 'SpecReg'
                [AllFramesFTrealign, MRS_struct] = Spectral_Registration(MRS_struct,0);
                MRS_struct.spec.AllFramesFTrealign = AllFramesFTrealign;
            case 'SpecRegDual'
                %Dual-channel Spectral Registration is applied separately to ON and OFF and they are coregistered after...
                [AllFramesFTrealign, MRS_struct] = Spectral_Registration(MRS_struct,0,1);
                MRS_struct.spec.AllFramesFTrealign = AllFramesFTrealign;
            case 'SpecRegHERMES' % MM (170703)
                [AllFramesFTrealign, MRS_struct] = Spectral_Registration_HERMES(MRS_struct);
                MRS_struct.spec.AllFramesFTrealign = AllFramesFTrealign;
            case 'none' % GO (180224)
                % do nothing
                MRS_struct.spec.AllFramesFTrealign = AllFramesFTrealign;
        end % end of switch for alignment
        
        % Separate ON/OFF data and generate DIFF spectra
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
            
            % Edit-OFF,-OFF spectrum (for Cr referencing) (MM: 180725)
            OFF_OFF = mean(AllFramesFTrealign(:,all([MRS_struct.fids.ON_OFF' MRS_struct.fids.ON_OFF2']==0,2) & MRS_struct.out.reject(:,ii)==0), 2);
            
            MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).off_off(ii,:)  = OFF_OFF;
            MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target2)).off_off(ii,:) = OFF_OFF;
            
            % Remove residual water from diff and diff_noalign spectra using HSVD -- GO & MGSaleh 2016
            if MRS_struct.p.water_removal
                
                fprintf('\nFiltering out residual water signal...\n');
                
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
        GannetPlotPrePostAlign(MRS_struct, vox, ii, kk);
        if MRS_struct.p.HERMES
            title({'Edited Spectra';'(pre- and post-alignment)'});
        else
            title({'Edited Spectrum';'(pre- and post-alignment)'});
        end
        xlabel('ppm');
        set(gca,'YTick',[]);
        
        % Top right
        hb = subplot(2,2,2);
        rejectframesplot = (1./MRS_struct.out.reject(:,ii).') .*  MRS_struct.spec.F0freq(ii,:);
        plot(1:size(FullData,2), MRS_struct.spec.F0freq(ii,:)', '-', 1:size(FullData,2), rejectframesplot, 'ro');
        set(gca,'XLim',[0 size(FullData,2)]);
        xlabel('average'); ylabel('\omega_0');
        if MRS_struct.p.HERMES
            title('Cr Frequency');
        else
            title('Water Frequency');
        end
        
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
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
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
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
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
            pdfname = fullfile('GannetLoad_output', [fullpath '_' vox{kk} '_load.pdf']); % MM (180112)
        else
            pdfname = fullfile('GannetLoad_output', [metabfile_nopath '_' vox{kk} '_load.pdf']); % MM (180112)
        end
        saveas(h, pdfname);
        
        
        % Export the processed data into an SDAT file
        if MRS_struct.p.sdat 
            if strcmpi(MRS_struct.p.vendor,'Philips')
                % Set up filenames
                sdat_G_name = ['GannetLoad_output/' metabfile_nopath  '_' vox{kk} '_G.sdat'];
                spar_G_name = ['GannetLoad_output/' metabfile_nopath  '_' vox{kk} '_G.spar'];
                % Make file copies for SDAT/SPAR files
                copyfile(metabfile{ii},sdat_G_name);
                sparname = [metabfile{ii}(1:end-4) MRS_struct.p.spar_string];
                copyfile(sparname,spar_G_name);
                % Write DIFF data into the SDAT file
                sdat_diff_out = conj(ifft(fftshift(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).diff(ii,:),2),[],2));
                sdat_diff_out = sdat_diff_out(1:MRS_struct.p.npoints(ii));
                % Also write out OFF data
                sdat_off_out = conj(ifft(fftshift(MRS_struct.spec.(vox{kk}).(sprintf('%s',MRS_struct.p.target)).off(ii,:),2),[],2));
                sdat_off_out = sdat_off_out(1:MRS_struct.p.npoints(ii));
                fileid  = fopen(sdat_G_name,'w','ieee-le');
                ff(:,1:2:2*MRS_struct.p.npoints(ii)) = real(sdat_diff_out);
                ff(:,2:2:2*MRS_struct.p.npoints(ii)) = imag(sdat_diff_out);
                gg(:,1:2:2*MRS_struct.p.npoints(ii)) = real(sdat_off_out);
                gg(:,2:2:2*MRS_struct.p.npoints(ii)) = imag(sdat_off_out);
                fwriteVAXD(fileid,[ff.' gg.'],'float');
                fclose(fileid);
            else
                warning('Only Philips SDAT files can be exported! No data exported.');
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
        
        % Save MRS_struct as mat file
        if ii == numscans && MRS_struct.p.mat
            % Set up filename
            mat_name = ['GannetLoad_output/MRS_struct_' vox{kk} '.mat'];
            save(mat_name,'MRS_struct');
        end
        
    end % end of output loop over voxels
    
end % end of load-and-processing loop over datasets

end



