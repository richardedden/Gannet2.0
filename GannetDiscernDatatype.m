function MRS_struct = GannetDiscernDatatype(filename,MRS_struct)
%Use the file ending to determine file type

lastchar=filename;
lastchar=lastchar((end-1):end);

if(strcmpi(lastchar,'.7'))
    MRS_struct.vendor = 'GE';
    MRS_struct.Reference_compound='H2O';
elseif(strcmpi(lastchar,'AT'))
    MRS_struct.vendor = 'Philips';
    if(strcmp(lastchar,'AT'))
       MRS_struct.spar_string='SPAR';
    else
        MRS_struct.spar_string='spar';
    end
elseif(strcmpi(lastchar,'TA'))
    MRS_struct.vendor = 'Philips_data';
elseif(strcmpi(lastchar,'DA'))
    MRS_struct.vendor = 'Siemens';
else
    error('Unrecognised filetype: should end .7 .SDAT or .RDA')
end
end
