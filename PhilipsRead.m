function [ MRS_struct ] = PhilipsRead(MRS_struct, fname, fname_water )
% RE/CJE Parse SPAR file for header info

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PhilipsRead is designed to handle 'off-first' data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 110825

   % work out data header name
   sparname = [fname(1:(end-4)) MRS_struct.p.spar_string]
   sparheader = textread(sparname, '%s');
   sparidx=find(ismember(sparheader, 'samples')==1);
   MRS_struct.p.npoints = str2num(sparheader{sparidx+2});
   sparidx=find(ismember(sparheader, 'rows')==1);
   MRS_struct.p.nrows = str2num(sparheader{sparidx+2});
   
   sparidx=find(ismember(sparheader, 'averages')==1);
   MRS_struct.p.Navg(MRS_struct.ii) = MRS_struct.p.nrows * str2num(sparheader{sparidx+2}); %raee 120228 fixing sdat NAvg
   %MRS_struct.p.Navg(MRS_struct.ii) = MRS_struct.p.nrows; %Trial SDAT might be average not sum.
   sparidx=find(ismember(sparheader, 'repetition_time')==1);
   MRS_struct.p.TR = str2num(sparheader{sparidx+2});
   
   sparidx=find(ismember(sparheader, 'sample_frequency')==1);
   MRS_struct.p.sw = str2num(sparheader{sparidx+2});
   
   %Pull in geometry information for GannetCoRegister
   %NP edit 150714; Below section is broken and now commented out as a temporary fix.
   %sparheadinfo is not set and may have to be sparheader instead. Also,
   %not looped yet (added edit0

   %for ii = 
% sparidx=find(ismember(sparheadinfo, 'ap_size')==1);
% MRS_struct.p.voxsize(ii,2) = str2num(sparheadinfo{sparidx+2});
% sparidx=find(ismember(sparheadinfo, 'lr_size')==1);
% MRS_struct.p.voxsize(ii,1) = str2num(sparheadinfo{sparidx+2});
% sparidx=find(ismember(sparheadinfo, 'cc_size')==1);
% MRS_struct.p.voxsize(ii,3) = str2num(sparheadinfo{sparidx+2});
% 
% sparidx=find(ismember(sparheadinfo, 'ap_off_center')==1);
% MRS_struct.p.voxoff(ii,2) = str2num(sparheadinfo{sparidx+2});
% sparidx=find(ismember(sparheadinfo, 'lr_off_center')==1);
% MRS_struct.p.voxoff(ii,1) = str2num(sparheadinfo{sparidx+2});
% sparidx=find(ismember(sparheadinfo, 'cc_off_center')==1);
% MRS_struct.p.voxoff(ii,3) = str2num(sparheadinfo{sparidx+2});
% 
% sparidx=find(ismember(sparheadinfo, 'ap_angulation')==1);
% MRS_struct.p.voxang(ii,2) = str2num(sparheadinfo{sparidx+2});
% sparidx=find(ismember(sparheadinfo, 'lr_angulation')==1);
% MRS_struct.p.voxang(ii,1) = str2num(sparheadinfo{sparidx+2});
% sparidx=find(ismember(sparheadinfo, 'cc_angulation')==1);
% MRS_struct.p.voxang(ii,3) = str2num(sparheadinfo{sparidx+2});

   
   
   MRS_struct.fids.data = SDATreadMEGA(fname, MRS_struct.p.npoints, MRS_struct.p.nrows);
   if nargin>2
       % work out data header name
       sparnameW = [fname_water(1:(end-4)) 'spar'];
       sparheaderW = textread(sparnameW, '%s');
       sparidxW=find(ismember(sparheaderW, 'averages')==1);
       %MRS_struct.p.Nwateravg = str2num(sparheaderW{sparidxW+2});
       MRS_struct.p.Nwateravg = 1; %SDAT IS average not sum.
   end
   %undo time series phase cycling to match GE
   corrph = ones(size(MRS_struct.fids.data));
   for jj=1:size(MRS_struct.fids.data,2)
    corrph(:,jj) = corrph(:,jj) * (-1).^jj;
   end
   
   MRS_struct.fids.data = MRS_struct.fids.data .* corrph;
   %Re-introduce initial phase step...
   MRS_struct.fids.data = MRS_struct.fids.data .*repmat(conj(MRS_struct.fids.data(1,:))./abs(MRS_struct.fids.data(1,:)),[MRS_struct.p.npoints 1]);
   %Philips data appear to be phased already (ideal case)
   
   %MRS_struct.fids.data = -conj(MRS_struct.fids.data); %RE 110728 - empirical factor to scale 'like GE'
   %If on-first - use the above...
   %If off-frst use:  FOR NOW, THIS DIDNT HELP...
   MRS_struct.fids.data = conj(MRS_struct.fids.data); %RE 110728 - empirical factor to scale 'like GE'
   if nargin>2
       % load water data
       MRS_struct.fids.data_water = SDATread(fname_water, MRS_struct.p.npoints);
       MRS_struct.fids.data_water(1,5);
       MRS_struct.fids.data_water = MRS_struct.fids.data_water.*conj(MRS_struct.fids.data_water(1))./abs(MRS_struct.fids.data_water(1));
       MRS_struct.out.phase_water = conj(MRS_struct.fids.data_water(1))./abs(MRS_struct.fids.data_water(1));
       %MRS_struct.fids.data = MRS_struct.fids.data.* MRS_struct.out.phase_water;
   end
end

