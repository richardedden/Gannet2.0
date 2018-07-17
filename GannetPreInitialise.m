function MRS_struct = GannetPreInitialise(MRS_struct)

% Some of these parameters will be parsed from the MRS data file header

% Acquisition Parameters
    MRS_struct.p.seqorig = 'JHU'; % origin of Philips patch; options are 'JHU' or 'Philips'
    MRS_struct.p.target = 'GABAGlx'; % signal(s) to fit; options are 'GABAGlx', 'GSH' or 'Lac'
    MRS_struct.p.target2 = 'GSH'; % applies to HERMES data only; options are 'GSH' or 'Lac'
    MRS_struct.p.ONOFForder = 'offfirst'; % order of editing pulses; options are 'onfirst' or 'offfirst'
    MRS_struct.p.Water_Positive = 1; % for Philips MOIST water suppression, set to 0
    
% Analysis Parameters
    MRS_struct.p.LB = 3; % line-broadening (in Hz)
    MRS_struct.p.water_phase_correction = 1; % perform phase correction; 1 = YES
    MRS_struct.p.data_phase_correction = 0; % perform phase correction; 1 = YES
    MRS_struct.p.water_removal = 1; % remove residual water using HSVD in HERMES (recommended for GABA/GSH editing) or GSH-edited MEGA-PRESS data; 1 = YES
    MRS_struct.p.AlignTo = 'SpecReg'; % options are 'SpecReg' (recommended for MEGA-PRESS), 'SpecRegHERMES' (recommended for HERMES and GSH-edited MEGA-PRESS), 'Cr', 'Cho', 'NAA', 'H2O' or 'CrOFF'
    MRS_struct.p.Vox = {'vox1','vox2'}; % for naming voxels acquired by PRIAM, e.g.: {'anterior','posterior'}, {'right','left'}, etc.
    MRS_struct.p.GSH_model = 'FiveGauss'; % choice of model for fitting GSH;
                                          % options are 'FiveGauss' (recommended for medium-TE HERMES) or 'SixGauss' (recommended for long-TE MEGA-PRESS)
    
% Flags
    MRS_struct.p.HERMES  = 0; % 1 = YES, 0 = NO (for MEGA-PRESS)
    MRS_struct.p.PRIAM   = 0; % 1 = YES, 0 = NO
    MRS_struct.p.phantom = 0; % 1 = YES, 0 = NO (for in vivo data)
    MRS_struct.p.mat     = 0; % 1 = YES, save MRS_struct as .mat file
    MRS_struct.p.sdat    = 0; % 1 = YES, save processed difference spectrum as .sdat file (only for Philips SDAT MEGA-PRESS datasets)
    MRS_struct.p.csv     = 0; % 1 = YES, extract useful data from MRS_struct and export to .csv file (applies to GannetFit, GannetSegment and GannetQuantify)
    
end