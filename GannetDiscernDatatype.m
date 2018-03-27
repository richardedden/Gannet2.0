function MRS_struct = GannetDiscernDatatype(filename,MRS_struct)
%Use the file ending to determine file type

lastchar=filename;

last2char=lastchar((end-1):end);
last3char=lastchar((end-2):end);
last4char=lastchar((end-3):end);

if strcmpi(last2char,'.7')
    MRS_struct.p.vendor = 'GE';
    MRS_struct.p.Reference_compound = 'H2O';
elseif strcmpi(last4char,'SDAT')
    MRS_struct.p.vendor = 'Philips';
    if strcmp(last4char,'SDAT')
        MRS_struct.p.spar_string = 'SPAR';
    else
        MRS_struct.p.spar_string = 'spar';
    end
elseif strcmpi(last2char,'TA')
    MRS_struct.p.vendor = 'Philips_data';
elseif(strcmpi(last4char,'.RAW')) % GO 11/01/2016
    MRS_struct.p.vendor = 'Philips_raw'; % GO 11/01/2016
elseif strcmpi(last3char,'RDA')
    MRS_struct.p.vendor = 'Siemens_rda';
elseif strcmpi(last4char,'.DAT')
    MRS_struct.p.vendor = 'Siemens_twix';
elseif(strcmpi(last4char,'.IMA')) % GO 11/11/2016
    MRS_struct.p.vendor = 'Siemens_dicom'; % GO 11/11/2016
elseif(strcmpi(last4char,'.DCM')) % GO 11/30/2016
    MRS_struct.p.vendor = 'dicom'; % GO 11/30/2016
else
    error('Unrecognised filetype: should end .7 .SDAT .DATA .RAW .RDA .IMA .DCM or .DAT'); % GO 11/11/2016
end

end
