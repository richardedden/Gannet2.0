function MRS_struct=GannetPreInitialise(MRS_struct)

% Some of these parameters will be overwritten by correct values are stored
% in the data headers.


%Acquisition Parameters
%    MRS_struct.sw=5000;  % sw taken from header for all formats except Philips .data
    MRS_struct.sw=2000; %This should be parsed from headers where possible    
    MRS_struct.TR=2000;%This should be parsed from headers where possible
    MRS_struct.TE=68; %This should be parsed from headers where possible
    MRS_struct.LarmorFreq=127.8; %This should be parsed from headers where possible
    MRS_struct.ONOFForder='offfirst';
    %Options are MRS_struct.ONOFForder='onfirst' or 'offfirst';
    MRS_struct.Water_Positive=1;
    
%Analysis Parameters
    MRS_struct.LB = 3;
    MRS_struct.ZeroFillTo = 32768;
    %AlignTo planned options: Cr; Cho; NAA; H20; CrOFF
    MRS_struct.AlignTo = 'SpecReg';



end
