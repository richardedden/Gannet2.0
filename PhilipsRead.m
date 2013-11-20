function [ MRS_struct ] = PhilipsRead(MRS_struct, fname, fname_water )
% RE/CJE Parse SPAR file for header info

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PhilipsRead is designed to handle 'off-first' data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 110825

   % work out data header name
   sparname = [fname(1:(end-4)) MRS_struct.spar_string]
   sparheader = textread(sparname, '%s');
   sparidx=find(ismember(sparheader, 'samples')==1);
   MRS_struct.npoints = str2num(sparheader{sparidx+2});
   sparidx=find(ismember(sparheader, 'rows')==1);
   MRS_struct.nrows = str2num(sparheader{sparidx+2});
   
   sparidx=find(ismember(sparheader, 'averages')==1);
   MRS_struct.Navg(MRS_struct.ii) = MRS_struct.nrows * str2num(sparheader{sparidx+2}); %raee 120228 fixing sdat NAvg
   %MRS_struct.Navg(MRS_struct.ii) = MRS_struct.nrows; %Trial SDAT might be average not sum.
   sparidx=find(ismember(sparheader, 'repetition_time')==1);
   MRS_struct.TR = str2num(sparheader{sparidx+2});
   
   sparidx=find(ismember(sparheader, 'sample_frequency')==1);
   MRS_struct.sw = str2num(sparheader{sparidx+2});
   
   MRS_struct.data = SDATreadMEGA(fname, MRS_struct.npoints, MRS_struct.nrows);
   if nargin>2
       % work out data header name
       sparnameW = [fname_water(1:(end-4)) 'spar'];
       sparheaderW = textread(sparnameW, '%s');
       sparidxW=find(ismember(sparheaderW, 'averages')==1);
       %MRS_struct.Nwateravg = str2num(sparheaderW{sparidxW+2});
       MRS_struct.Nwateravg = 1; %SDAT IS average not sum.
   end
   %undo time series phase cycling to match GE
   corrph = ones(size(MRS_struct.data));
   for jj=1:size(MRS_struct.data,2)
    corrph(:,jj) = corrph(:,jj) * (-1).^jj;
   end
   
   MRS_struct.data = MRS_struct.data .* corrph;
   %Re-introduce initial phase step...
   MRS_struct.data = MRS_struct.data .*repmat(conj(MRS_struct.data(1,:))./abs(MRS_struct.data(1,:)),[MRS_struct.npoints 1]);
   %Philips data appear to be phased already (ideal case)
   
   %MRS_struct.data = -conj(MRS_struct.data); %RE 110728 - empirical factor to scale 'like GE'
   %If on-first - use the above...
   %If off-frst use:  FOR NOW, THIS DIDNT HELP...
   MRS_struct.data = conj(MRS_struct.data); %RE 110728 - empirical factor to scale 'like GE'
   if nargin>2
       % load water data
       MRS_struct.data_water = SDATread(fname_water, MRS_struct.npoints);
       MRS_struct.data_water(1,5)
       MRS_struct.data_water = MRS_struct.data_water.*conj(MRS_struct.data_water(1))./abs(MRS_struct.data_water(1));
       MRS_struct.phase_water = conj(MRS_struct.data_water(1))./abs(MRS_struct.data_water(1));
       %MRS_struct.data = MRS_struct.data.* MRS_struct.phase_water;
   end
end

