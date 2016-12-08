function GEDeIdentify(fnames)
% GEDeIdentify
%   Reads GE P-files and removes participant/exam information from the
%   header. Newly de-identified P-files are then output, with filenames
%   appended with '_noID'. The original files are not overwritten.
%
%   If GEDeIdentify has already been run, and de-identified files are found
%   within the directory, the user will be asked whether these files should
%   be overwritten.
%
%   NOTE: The user must make sure that filenames themselves do not contain
%   information that can personally identify participants. This function
%   will only de-identify the content of the P-file header.
%
%   Usage:
%       GEDeIdentify, by itself, de-identifies all P-files found within the
%       current directory.
%
%       GEDeIdentify(fnames) de-identifies P-files listed in the cell array
%       fnames.
%
%   Example:
%       c = {'P15221.7', 'P17981.7'};
%       GEDeIdentify(c);

%   Author: Mark Mikkelsen (Johns Hopkins University, Aug. 2016)
%
%   Acknowledgements:
%       This code is heavily based on plotp.m and pfile_anon.c written by
%       Fred Frigo (Marquette University). The author is also grateful to
%       Ralph Noeske (GE Healthcare) for providing details regarding the
%       P-file header.
%
%   Version history:
%       2016-08-09: + Function created


if nargin < 1 % De-identify all P-files in current directory
    
    flist = dir('*.7');
    for ii = 1:length(flist)
        fnames(ii) = cellstr(flist(ii).name);
    end
    
    nArgs = nargin;
    [exitFunc, fnames] = CheckForOutput(nArgs, fnames);
    if exitFunc
        return
    end
    
else % De-identify P-files user has listed in fnames
    
    % Check if filenames include a .7 extension
    for ii = 1:length(fnames)
        ext = fnames{ii}(end-1:end);
        assert(strcmpi(ext, '.7'), ['The filename ' fnames{ii} ' does not include a .7 extension.']);
    end
    
    % Check if files can be found
    for ii = 1:length(fnames)
        assert(any(exist(fnames{ii}, 'file')), ...
            ['The file ' fnames{ii} ' cannot be found.' ...
            ' Check spelling of filenames (P-files must include an extension in their filename).' ...
            ' Also check that you are in the right directory.']);
    end
    
    nArgs = nargin;
    [exitFunc, fnames] = CheckForOutput(nArgs, fnames);
    if exitFunc
        return
    end
    
end

% Read P-files and remove participant/exam information; save new P-files
for ii = 1:length(fnames)
    
    pfile_fid = fopen(fnames{ii}, 'r', 'ieee-be');
    pfile_fid_noID = fopen([fnames{ii}(1:end-2) '_noID' fnames{ii}(end-1:end)], 'w', 'ieee-be');
    
    % Check if P-file can be found
    if pfile_fid == -1
        fclose(pfile_fid);
        fclose(pfile_fid_noID);
        error(['The file ' fnames{ii} ' cannot be found.' ...
            ' Check spelling of filenames (P-files must include an extension in their filename).' ...
            ' Also check that you are in the right directory.']);
    end
    
    % Determine P-file version
    [pfile_fid, pfile_fid_noID, pfile_hdr_sz, pfile_exam_offset, pfile_series_offset] = ...
        VersionCheck(fnames{ii}, pfile_fid, pfile_fid_noID);
    
    % Determine how big the P-file is
    file_info = dir(fnames{ii});
    pfile_sz = file_info.bytes;
    
    % Read the entire P-file
    frewind(pfile_fid);
    pfile = fread(pfile_fid, pfile_sz, 'integer*2');
    
    % Read some header information
    frewind(pfile_fid);
    hdr_val = fread(pfile_fid, 102, 'integer*2');
    n_echoes = hdr_val(36);
    point_sz = hdr_val(42);
    n_points = hdr_val(52);
    n_rows = hdr_val(53);
    n_receivers = (hdr_val(102) - hdr_val(101)) + 1;
            
    % Spectro prescan P-files
    if n_points == 1 && n_rows == 1
        n_points = 2048;
    end
    
    % Compute size (in bytes) of raw data
    data_elements = n_points * 2;
    my_frame = 1;    
    total_frames = (n_rows - my_frame + 1) * n_echoes;
    data_elements = data_elements * total_frames * n_receivers;
    
    % Read the raw data (using the appropriate precision)
    fseek(pfile_fid, pfile_hdr_sz, 'bof');
    if point_sz == 2
        pfile_raw_data = fread(pfile_fid, data_elements, 'integer*2');
    else
        pfile_raw_data = fread(pfile_fid, data_elements, 'integer*4');
    end    
    
    % Write all P-file bytes into new P-file
    frewind(pfile_fid_noID);
    fwrite(pfile_fid_noID, pfile, 'integer*2');
    
    % Overwrite EXAMDATATYPE structure that contains participant/exam
    % information with zeros in new P-file
    exam_data_type_sz = pfile_series_offset - pfile_exam_offset;
    fseek(pfile_fid_noID, pfile_exam_offset, 'bof');
    fwrite(pfile_fid_noID, zeros(exam_data_type_sz,1), 'integer*2');
    
    % Write raw data to new P-file (using the appropriate precision)
    fseek(pfile_fid_noID, pfile_hdr_sz, 'bof');
    if point_sz == 2
        fwrite(pfile_fid_noID, pfile_raw_data, 'integer*2');
    else
        fwrite(pfile_fid_noID, pfile_raw_data, 'integer*4');
    end
    
    fclose(pfile_fid);
    fclose(pfile_fid_noID);
    
end


function [exitFunc, fnames] = CheckForOutput(nArgs, fnames)
% Check if any de-identified files have already been output and ask user if
% they want to overwrite them

exitFunc = 0;

if nArgs < 1
    
    if any(~cellfun('isempty', strfind(fnames, '_noID')))
        resp = input('\nDe-identified files found in the directory! Proceed and overwrite? [y/n]: ','s');
        if strcmpi(resp, 'y')
            disp('Overwriting...');
        elseif strcmpi(resp, 'n')
            disp('Exiting...');
            exitFunc = 1;
            return
        end
        
        ind1 = strfind(fnames, '_noID');
        ind2 = find(~cellfun('isempty', ind1));
        fnames(ind2) = []; %#ok<FNDSB>
    end
    
else
    
    for ii = 1:length(fnames)
        fnames_noID{ii} = [fnames{ii}(1:end-2) '_noID' fnames{ii}(end-1:end)]; %#ok<AGROW>
    end
    
    if any(cellfun(@exist, fnames_noID))
        resp = input('\nDe-identified files found in the directory! Proceed and overwrite? [y/n]: ','s');
        if strcmpi(resp, 'y')
            disp('Overwriting...');
        elseif strcmpi(resp, 'n')
            disp('Exiting...');
            exitFunc = 1;
            return
        end
        
    end
    
end


function [pfile_fid, pfile_fid_noID, pfile_hdr_sz, pfile_exam_offset, pfile_series_offset] = ...
    VersionCheck(fnames, pfile_fid, pfile_fid_noID)
% Determine P-file version

frewind(pfile_fid);
rdbm_rev_num = fread(pfile_fid, 1, 'real*4');

if rdbm_rev_num == 7.0 % Signa 8.X (LX)
    pfile_hdr_sz = 39984;
    pfile_exam_offset = pfile_hdr_sz - 1040 - 1028 - 1044;
    pfile_series_offset = pfile_hdr_sz - 1028 - 1044;
elseif rdbm_rev_num == 8.0 % Cardiac / MGD
%     pfile_hdr_sz = 60464;
    fclose(pfile_fid);
    fclose(pfile_fid_noID);
    error('GEDeIdentify does not yet support P-file header revision: %f', rdbm_rev_num);
elseif rdbm_rev_num > 5.0 && rdbm_rev_num < 6.0 % Signa 5.5
%     pfile_hdr_sz = 39940;
    fclose(pfile_fid);
    fclose(pfile_fid_noID);
    error('GEDeIdentify does not yet support P-file header revision: %f', rdbm_rev_num);
else
    % In 11.0 and later, the header and data are stored as little-endian
    fclose(pfile_fid);
    fclose(pfile_fid_noID);
    pfile_fid = fopen(fnames, 'r', 'ieee-le');
    pfile_fid_noID = fopen([fnames(1:end-2) '_noID' fnames(end-1:end)], 'w', 'ieee-le');
    
    frewind(pfile_fid);
    rdbm_rev_num = fread(pfile_fid, 1, 'real*4');
    
    if rdbm_rev_num == 9.0  % Excite2 11.0
        pfile_hdr_sz = 61464;
        pfile_exam_offset = pfile_hdr_sz - 1040 - 1536 - 1536;
        pfile_series_offset = pfile_hdr_sz - 1536 - 1536;
    elseif rdbm_rev_num >= 11.0 && rdbm_rev_num < 100.0  % 12.0 and later
        fseek(pfile_fid, 1468, 'bof');
        pfile_hdr_sz = fread(pfile_fid, 1, 'integer*4');
        fseek(pfile_fid, 1496, 'bof');
        pfile_exam_offset = fread(pfile_fid, 1, 'integer*4');
        fseek(pfile_fid, 1500, 'bof');
        pfile_series_offset = fread(pfile_fid, 1, 'integer*4');
    else
        fclose(pfile_fid);
        fclose(pfile_fid_noID);
        error('Invalid P-file header revision: %f', rdbm_rev_num);
    end
    
end


