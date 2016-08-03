function PhilipsDeIdentify(fnames)
% PhilipsDeIdentify(fnames)
%   Reads Philips SPAR files and removes participant information. New
%   de-identified SPAR files are then output, with filenames appended with
%   '_noID'. The original files are not overwritten.
%
%   NOTE: The user must make sure that filenames themselves do not contain
%   information that can personally identify participants. This function
%   will only de-identify the content of the SPAR file.
%
%   Usage:
%     fnames = Cell array containing SPAR filenames
%
%   Example:
%     PhilipsDeIdentify({'S01_gaba_7_2_raw_act.SPAR', 'S01_gaba_7_2_raw_ref.SPAR'});
%
%   Author: Dr. Mark Mikkelsen (Johns Hopkins University, 2016-08-03)


% Check if filenames include a .SPAR/.spar extension
for ii = 1:length(fnames)
    ext = fnames{ii}(end-4:end);
    assert(strcmpi(ext,'.spar'), ...
        ['The filename ' fnames{ii} ' (' num2str(ii) ')' ' in ' inputname(1) ...
        ' does not include a .SPAR/.spar extension.']);
end

% Check if files can be found
for ii = 1:length(fnames)
    assert(any(exist(fnames{ii},'file')), ...
        ['The file ' fnames{ii} ' (' num2str(ii) ')' ' cannot be found.' ...
        ' Check spelling of filenames in ' inputname(1) ...
        ' (SPAR files must include an extension in their filename).' ...
        ' Also check that you are in the right directory.']);
end

% Read SPAR files and remove participant information; save new SPAR files
for ii = 1:length(fnames)
    
    spar_fid = fopen(fnames{ii}, 'r');
    spar_fid_noID = fopen([fnames{ii}(1:end-5) '_noID' fnames{ii}(end-4:end)], 'w');
    
    tline = fgetl(spar_fid);
    while ischar(tline)
        if any(strfind(tline, 'examination_name'));
            tline = 'examination_name : ';
        elseif any(strfind(tline, 'patient_name'));
            tline = 'patient_name : ';
        elseif any(strfind(tline, 'patient_birth_date'));
            tline = 'patient_birth_date : ';
        end
        fprintf(spar_fid_noID, '%s\n', tline);
        tline = fgetl(spar_fid);
    end
    
    fclose(spar_fid);
    fclose(spar_fid_noID);
    
end



