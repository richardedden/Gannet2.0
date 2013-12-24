function MRS_struct=GannetPreInitialise(MRS_struct)

% Some of these parameters will be overwritten by correct values are stored
% in the data headers.


%Acquisition Parameters
%    MRS_struct.sw=5000;  % sw taken from header for all formats except Philips .data
    MRS_struct.sw=1600; %This should be parsed from headers where possible
    MRS_struct.npoints=2048; %This is twice the acquired points for TWIX data;
    %This should be parsed from headers where possible
    MRS_struct.TR=2000;%This should be parsed from headers where possible
    MRS_struct.TE=68; %This should be parsed from headers where possible
    MRS_struct.LarmorFreq=123.0; %This should be parsed from headers where possible
    %In general, LarmorFreq is 127.8 on Philips,
    MRS_struct.target='GABA'; %Other option is GSH
    MRS_struct.ONOFForder='onfirst';
    %Options are MRS_struct.ONOFForder='onfirst' or 'offfirst';
    MRS_struct.Water_Positive=1;
    
%Analysis Parameters
    MRS_struct.LB = 3;
    MRS_struct.ZeroFillTo = 32768;
    %AlignTo planned options: Cr; Cho; NAA; H20; CrOFF
    MRS_struct.AlignTo = 'SpecReg';
end
