function MRS_struct = GannetPreInitialise(MRS_struct)

% Some of these parameters will be parsed from the data file header

% Acquisition Parameters
    MRS_struct.p.sw = []; % parsed from header
    MRS_struct.p.npoints = []; % parsed from header
    MRS_struct.p.TR = []; % parsed from header
    MRS_struct.p.TE = []; % parsed from header
    MRS_struct.p.LarmorFreq = []; % parsed from header
    MRS_struct.p.Nwateravg = []; % parsed from header
    MRS_struct.p.target = 'GABAGlx'; % for either MEGA-PRESS or HERMES; options are 'GABAGlx' or 'GSH'
    MRS_struct.p.target2 = 'GSH'; % for HERMES only; options are 'GSH' or 'Lac'
    MRS_struct.p.ONOFForder = 'offfirst'; % options are 'onfirst' or 'offfirst'
    MRS_struct.p.Water_Positive = 1; % for Philips MOIST ws, set to 0
    % Siemens header information differs between versions; switch for different versions
    MRS_struct.p.Siemens_type = 1; % 1 = TIM TRIO WIP; 2 = Near seq; 3 = Skyra WIP; 4 = Prisma (VD13C); 5 = Prisma (Minnesota); 6 = Jamie's VE11B (Jena)
        
% Analysis Parameters
    MRS_struct.p.LB = 3; % line-broadening (in Hz)
    MRS_struct.p.water_phase_correction = 1; % perform phase correction; 1 = YES
    MRS_struct.p.data_phase_correction = 0; % perform phase correction; 1 = YES
    MRS_struct.p.water_removal = 1; % remove residual water in HERMES data using HSVD; 1 = YES
    MRS_struct.p.AlignTo = 'SpecReg'; % options are 'SpecReg' (recommended for MEGA-PRESS), 'SpecRegHERMES' (recommended for HERMES), 'Cr', 'Cho', 'NAA', 'H2O', 'CrOFF'
    MRS_struct.p.Vox = {'vox1','vox2'}; % for naming voxels acquired by PRIAM, e.g: 'anterior' and 'posterior', 'right' and 'left', etc.
    
% Flags
    MRS_struct.p.HERMES = 0; % 1 = YES, 0 = NO (MEGA-PRESS)
    MRS_struct.p.PRIAM  = 0; % 1 = YES, 0 = NO
    MRS_struct.p.mat    = 0; % 1 = YES, save MRS_struct as .mat file
    MRS_struct.p.sdat   = 0; % 1 = YES, save MRS_struct as .sdat file
    
end