function MRS_struct=GannetPreInitialise(MRS_struct)

% Some of these parameters will be overwritten by correct values are stored
% in the data headers.
%mgs

%Acquisition Parameters
%    MRS_struct.sw=5000;  % sw taken from header for all formats except Philips .data
   % MRS_struct.p.sw=5000; %This should be parsed from headers where possible
%    MRS_struct.p.npoints=4096; %This is twice the acquired points for TWIX data;
    

    %This should be parsed from headers where possible
    MRS_struct.p.TR=2000;%This should be parsed from headers where possible
    MRS_struct.p.TE=68; %This should be parsed from headers where possible. It is parsed from the header -- MGSaleh 2016.
    MRS_struct.p.LarmorFreq=127; %This should be parsed from headers where possible. It is parsed from the header -- MGSaleh 2016.
    %In general, LarmorFreq is 127.8 on Philips,
    MRS_struct.p.Nwateravg = 8; %Needed for GE
    
    
    MRS_struct.p.target='GABAGlx'; % Options for both MEGA-PRESS or HERMES (target 1). Options are only 'GABAGlx' or 'GSH' -- MGSaleh 2016
    MRS_struct.p.target2='GSH';    % Options are 'Lac' or 'GSH' (if not specified in target 1) -- MGSaleh 2016
    
    %%(now also implemented) 'GABAGlx', 'GSH', 'Lac' (Lactate) in dev
    
    MRS_struct.p.ONOFForder='offfirst';
    %Options are MRS_struct.ONOFForder='onfirst' or 'offfirst';
    MRS_struct.p.Water_Positive=1; %For Philips MOIST ws, set to 0.
    %Siemens header information differs between versions
    %switch for different versions
    MRS_struct.p.Siemens_type = 1; %1 = TIM TRIO WIP 2 = Near seq 3 =Skyra WIP; 4=Prisma (VD13C); 5=Prisma(Minnesota)
    
    % A choice to perform phase correction on the water or not. 
    % Default = 1. Yes, perform the correction -- MGSaleh 2016    
    MRS_struct.p.water_phase_correction = 1;
    MRS_struct.p.data_phase_correction = 0;
    
    %Removing water using HSLVD -- GO and MGSaleh 2016
    MRS_struct.p.water_removal = 01;
    
    
%Analysis Parameters
    MRS_struct.p.LB = 3;
    MRS_struct.p.ZeroFillTo = 32768;
    %AlignTo planned options: Cr; Cho; NAA; H20; CrOFF
    MRS_struct.p.AlignTo = 'Cho'; %SpecReg default and recommended
    
%Flags
    MRS_struct.p.HERMES = 1;                % 1 = YES,  0 = NO (means MEGA-PRESS);  % Added by MGSaleh 2016 
    MRS_struct.p.PRIAM  = 0;                 % 1 = YES (mean dual voxels),  0 = NO;  % Added by MGSaleh 2016 
    MRS_struct.p.mat    = 0;                 % 1 = YES, save MRS_struct as .mat file 2016
    MRS_struct.p.sdat   = 0;                % 1 = YES, save MRS_struct as .sdat file 2016
end
