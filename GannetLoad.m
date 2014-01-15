function MRS_struct=GannetLoad(gabafile, waterfile)
%Gannet 2.0 GannetLoad
%Started by RAEE Nov 5, 2012

%Aim to make the GannetLoad more modular and easier to understand/edit, and
%especially to integrate the workflow for different filetypes more.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Work flow Summary
%   1.Pre-initialise
%   2. Determine data parameters from headers
%   3. Some Housekeeping
%   4. Load Data from files
%   5. Apply appropriate pre-processing
%   6. Output processed spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   0. Check the file list for typos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
missing=0;
for filecheck=1:length(gabafile)
    if(~exist(gabafile{filecheck}))
        disp(['The file ' gabafile{filecheck} ' (' num2str(filecheck) ')' ' is missing. Typo?'])
        missing=1;
    end
end
if(nargin > 1)
    for filecheck=1:length(waterfile)        
        if(~exist(waterfile{filecheck}))
            disp(['The file ' waterfile(filecheck) ' is missing. Typo?'])
            missing=1;
        end
    end
end
if missing
        error('Not all the files are there, so I give up.');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Pre-initialise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MRS_struct.versionload = 'G2 131202 raee';
MRS_struct.ii=0;
MRS_struct.gabafile=gabafile;
MRS_struct=GannetPreInitialise(MRS_struct);
%Check whether water data or not
if(nargin > 1)
    MRS_struct.waterfile = waterfile;
    MRS_struct.Reference_compound='H2O';
else
    MRS_struct.Reference_compound='Cr';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2. Determine data parameters from header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if iscell(gabafile) == 1 % it's a cell array, so work out the number of elements
        numpfiles=numel(gabafile);
        pfiles=gabafile;
    else
        numpfiles=1;  % it's just one pfile
        pfiles{1}=gabafile;
    end
    
MRS_struct=GannetDiscernDatatype(pfiles{1},MRS_struct);

    if(strcmpi(MRS_struct.vendor,'Siemens'))
        numpfiles = numpfiles/2;
    end

%%%%%%%%%%%%%%%%%%%%%%%%    
%   3. Some Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%

    % create dir for output
    if(exist('./MRSload_output','dir') ~= 7)
        mkdir MRSload_output
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%   4. Load Data from files
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:numpfiles    %Loop over all files in the batch (from gabafile)
    MRS_struct.ii=ii;
    
    switch MRS_struct.vendor
        case 'GE'
            Water_Positive=1;           %CHECK
            AlignTo = 2;           %CHECK
            MRS_struct = GERead(MRS_struct, gabafile{ii});
            da_xres = MRS_struct.npoints;
            da_yres = MRS_struct.nrows;
            WaterData = MRS_struct.data_water;
            MRS_struct.data = MRS_struct.data*MRS_struct.nrows/MRS_struct.Navg(ii);%I think GE does sum over NEX
            FullData = MRS_struct.data;
            ComWater = mean(WaterData,2);
            %Set up vector of which rows of .data are ONs and OFFs.
            switch MRS_struct.ONOFForder
                case 'onfirst'
                    MRS_struct.ON_OFF=repmat([1 0],[1 size(MRS_struct.data,2)/2]);
                case 'offfirst'
                    MRS_struct.ON_OFF=repmat([0 1],[1 size(MRS_struct.data,2)/2]);
            end
            totalframes = MRS_struct.nrows;
        case 'Siemens_twix'
            Water_Positive=1;           %CHECK
            if(exist('waterfile'))
                MRS_struct = SiemensTwixRead(MRS_struct, gabafile{ii}, waterfile{ii});
                WaterData = MRS_struct.data_water;
                ComWater = mean(WaterData,2);
            else
                MRS_struct = SiemensTwixRead(MRS_struct, gabafile{ii});
            end
            da_xres = MRS_struct.npoints;
            da_yres = MRS_struct.nrows;
            FullData = MRS_struct.data;
            %Set up vector of which rows of .data are ONs and OFFs.
            switch MRS_struct.ONOFForder
                case 'onfirst'
                    MRS_struct.ON_OFF=repmat([1 0],[1 size(MRS_struct.data,2)/2]);
                case 'offfirst'
                    MRS_struct.ON_OFF=repmat([0 1],[1 size(MRS_struct.data,2)/2]);
            end
            totalframes = MRS_struct.nrows;
        case 'Siemens'
            if(exist('waterfile'))    
                MRS_struct.Reference_compound='H2O';
                switch MRS_struct.ONOFForder
                    case 'offfirst'
                        MRS_struct = SiemensRead_RE(MRS_struct, gabafile{ii*2-1},gabafile{ii*2}, waterfile{ii});
                    case 'onfirst'
                        MRS_struct = SiemensRead_RE(MRS_struct, gabafile{ii*2},gabafile{ii*2-1}, waterfile{ii});
                end    
                MRS_struct.Nwateravg = 1;
                MRS_struct.phase{ii} = 0;
                MRS_struct.phase_firstorder(ii) = 0;
            else
                 MRS_struct.Reference_compound='Cr';
                switch MRS_struct.ONOFForder
                    case 'offfirst'
                        MRS_struct = SiemensRead_RE(MRS_struct, gabafile{ii*2-1},gabafile{ii*2});
                    case 'onfirst'
                        MRS_struct = SiemensRead_RE(MRS_struct, gabafile{ii*2},gabafile{ii*2-1});
                end    
             end
            da_xres = MRS_struct.npoints;
            da_yres = 1;
            totalframes = 1;
            FullData = MRS_struct.data;
            if(strcmp(MRS_struct.Reference_compound,'H2O'))
                WaterData = MRS_struct.data_water;
            end
            MRS_struct.LarmorFreq;
            % work out frequency scale 121106 (remving CSize)
            freqrange=MRS_struct.sw/MRS_struct.LarmorFreq;
            MRS_struct.freq=(MRS_struct.ZeroFillTo+1-(1:1:MRS_struct.ZeroFillTo))/MRS_struct.ZeroFillTo*freqrange+4.7-freqrange/2.0;
            MRS_struct.FreqPhaseAlign=0;
            %Data are always read in OFF then ON
            
            totalframes = 2;
            switch MRS_struct.ONOFForder
                case 'onfirst'
                    MRS_struct.ON_OFF=[1 0];
                    MRS_struct.ON_OFF=MRS_struct.ON_OFF(:).';%re
                case 'offfirst'
                    MRS_struct.ON_OFF=[0 1];
                    MRS_struct.ON_OFF=MRS_struct.ON_OFF(:).';
            end           
        case 'Philips'           
            if(exist('waterfile'))
                MRS_struct.Reference_compound='H2O';
            else
                 MRS_struct.Reference_compound='Cr';
            end
            %Need to set Water_Positive based on water signal
            if strcmpi(MRS_struct.Reference_compound,'H2O')
                MRS_struct = PhilipsRead(MRS_struct, gabafile{ii}, waterfile{ii});
                WaterData = MRS_struct.data_water;
            else
                MRS_struct = PhilipsRead(MRS_struct, gabafile{ii});
            end
            if MRS_struct.Water_Positive==0
                MRS_struct.data=-MRS_struct.data;
            end          
            da_xres = MRS_struct.npoints;
            da_yres = MRS_struct.nrows;
            totalframes = MRS_struct.nrows;
            FullData = MRS_struct.data;
            AlignTo = 2;           %CHECK
            switch MRS_struct.ONOFForder
                case 'onfirst'
                    MRS_struct.ON_OFF=repmat([1 0],[1 size(MRS_struct.data,2)/2]);
                case 'offfirst'
                    MRS_struct.ON_OFF=repmat([0 1],[1 size(MRS_struct.data,2)/2]);
            end
        case 'Philips_data'
            if(exist('waterfile'))    
                MRS_struct.Reference_compound='H2O';
                MRS_struct = PhilipsRead_data(MRS_struct, gabafile{ii},waterfile{ii});
            else
                 MRS_struct.Reference_compound='Cr';
                 MRS_struct = PhilipsRead_data(MRS_struct, gabafile{ii});
            end
            Water_Positive=1;           %CHECK
            if strcmpi(MRS_struct.Reference_compound,'H2O')
                WaterData = MRS_struct.data_water;
            end
            da_xres = MRS_struct.npoints;
            da_yres = MRS_struct.nrows*MRS_struct.Navg(ii);
            totalframes = MRS_struct.Navg(ii);
            FullData = MRS_struct.data;
            AlignTo = 1;           %CHECK
            switch MRS_struct.ONOFForder
                case 'onfirst'
                    MRS_struct.ON_OFF=repmat([1 0],[MRS_struct.Navg(ii)/MRS_struct.nrows MRS_struct.nrows/2]);
                    MRS_struct.ON_OFF=MRS_struct.ON_OFF(:).';
                case 'offfirst'
                    MRS_struct.ON_OFF=repmat([0 1],[MRS_struct.Navg(ii)/MRS_struct.nrows MRS_struct.nrows/2]);
                    MRS_struct.ON_OFF=MRS_struct.ON_OFF(:).';
            end
    end    %End of vendor switch loop for data load


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   5. Apply appropriate pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %There are some decisions to be made on what processing is applied to
    %what data

    %First steps
    MRS_struct.zf=MRS_struct.ZeroFillTo/MRS_struct.npoints;
    time=(1:1:size(FullData,1))/MRS_struct.sw;
    time_zeropad=(1:1:MRS_struct.ZeroFillTo)/(MRS_struct.sw);
    DataSize = size(FullData,2);

        % Finish processing water data. 
        if(strcmpi(MRS_struct.Reference_compound,'H2O'))     

            if(strcmpi(MRS_struct.vendor,'GE'))           %CHECK
                ComWater = mean(WaterData,2);           %CHECK
            elseif(strcmpi(MRS_struct.vendor,'Siemens'))       %CHECK
                ComWater = WaterData;           %CHECK
            elseif (strcmpi(MRS_struct.vendor,'Siemens_twix')) %CHECK
                ComWater = WaterData;
            else
                ComWater = WaterData.';           %CHECK
            end           %CHECK
            
            ComWater = ComWater.*exp(-(time')*MRS_struct.LB*pi);
            MRS_struct.spec.water(ii,:)=fftshift(fft(ComWater,MRS_struct.ZeroFillTo,1))';
        end %End of H20 reference loop

            FullData = FullData.* repmat( (exp(-(time')*MRS_struct.LB*pi)), [1 totalframes]);
            AllFramesFT=fftshift(fft(FullData,MRS_struct.ZeroFillTo,1),1);
            % work out frequency scale
            freqrange=MRS_struct.sw/MRS_struct.LarmorFreq;
            MRS_struct.freq=(MRS_struct.ZeroFillTo+1-(1:1:MRS_struct.ZeroFillTo))/MRS_struct.ZeroFillTo*freqrange+4.7-freqrange/2.0;

            %  Frame-by-frame Determination of max Frequency in spectrum (assumed water) maximum
            % find peak location for frequency realignment
            [FrameMax, FrameMaxPos] = max(AllFramesFT, [], 1);
            %Not always true that water starts at 4.68, if drift is rapid...
            water_off=abs(MRS_struct.freq-4.68);
            water_index=find(min(water_off)==water_off);
            % Determine Frame shifts
            FrameShift = FrameMaxPos - water_index;
            %Apply for Philips data, not for GE data (Why?)
            switch MRS_struct.vendor
                case 'GE'           %CHECK
                    AllFramesFTrealign=AllFramesFT;
                case {'Philips','Philips_data'}           %CHECK
                    for(jj=1:size(AllFramesFT,2))
                        AllFramesFTrealign(:,jj)=circshift(AllFramesFT(:,jj), -FrameShift(jj));             %CHECK - is this used????
                    end
                    %This quite possibly doesn't carry through, as it seems
                    %that the later stuff all starts with AllFramesFT, no
                    %AllFramesFTrealign.
            end %end of switch for Water max alignment p[re-initialisation

            MRS_struct.waterfreq(ii,:) = MRS_struct.freq(FrameMaxPos);%to be used for the output figure

            %Frame-by-Frame alignment
            switch MRS_struct.AlignTo
               case 'Cr'
                    [AllFramesFTrealign MRS_struct]=AlignUsingPeak(AllFramesFTrealign,MRS_struct);
               case 'Cr'
                    %AllFramesFTrealign=AlignUsingCr(AllFramesFTrealign,MRS_struct.ONOFForder,n);     
               case 'Cho'
                    %AllFramesFTrealign=AlignUsingCho(AllFramesFTrealign);
               case 'H20'
                   %AllFramesFTrealign=AlignUsingH2O(AllFramesFTrealign);
               case 'NAA'
                   [AllFramesFTrealign MRS_struct]=AlignUsingPeak(AllFramesFTrealign,MRS_struct);
               case 'SpecReg'
                    [AllFramesFTrealign MRS_struct] = Spectral_Registration(MRS_struct,0);
               case 'SpecRegDual'
                   %Dual-channel Spectral Registration is applied separately to ON and OFF and they are coregistered after... 
                   [AllFramesFTrealign MRS_struct] = Spectral_Registration(MRS_struct,0,1);
               end %end of switch for alignment target   

        
        %Separate ON/OFF data and generate SUM/DIFF (averaged) spectra.
        %In Gannet 2.0 Odds and Evens are explicitly replaced by ON and OFF        
        
        MRS_struct.spec.off(ii,:)=mean(AllFramesFTrealign(:,((MRS_struct.ON_OFF==0)'&(MRS_struct.reject(:,ii)==0))),2);
        MRS_struct.spec.on(ii,:)=mean(AllFramesFTrealign(:,((MRS_struct.ON_OFF==1)'&(MRS_struct.reject(:,ii)==0))),2);
        MRS_struct.spec.diff(ii,:)=(MRS_struct.spec.on(ii,:)-MRS_struct.spec.off(ii,:))/2; %Not sure whether we want a two here.
        MRS_struct.spec.diff_noalign(ii,:)=(mean(AllFramesFT(:,(MRS_struct.ON_OFF==1)),2)-mean(AllFramesFT(:,(MRS_struct.ON_OFF==0)),2))/2; %Not sure whether we want a two here.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   6. Build Gannet Output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        if(ishandle(101))
            close(101)
        end
        h=figure(101);
        set(h, 'Position', [100, 100, 1000, 707]);
        set(h,'Color',[1 1 1]);
        figTitle = ['GannetLoad Output'];
        set(gcf,'Name',figTitle,'Tag',figTitle, 'NumberTitle','off');
              
            %Top Left
            ha=subplot(2,2,1);
            Gannetplotprepostalign(MRS_struct,ii)
            x=title({'Edited Spectrum';'(pre- and post-align)'});
            set(gca,'YTick',[]);
            %Top Right
            hb=subplot(2,2,2);
            rejectframesplot = (1./MRS_struct.reject(:,ii).') .*  MRS_struct.waterfreq(ii,:);
            plot([1:DataSize], MRS_struct.waterfreq(ii,:)', '-', [1:DataSize], rejectframesplot, 'ro');
            set(gca,'XLim',[0 DataSize]);
            xlabel('time'); ylabel('\omega_0');
            title('Water Frequency, ppm');
             %Bottom Left
             hc=subplot(2,2,3);
             if strcmp(MRS_struct.AlignTo,'no')~=1
                CrFitLimLow=2.72;
                CrFitLimHigh=3.12;
                z=abs(MRS_struct.freq-CrFitLimHigh);
                lb=find(min(z)==z);
                z=abs(MRS_struct.freq-CrFitLimLow);
                ub=find(min(z)==z);
                CrFitRange=ub-lb;
                plotrealign=[ real(AllFramesFT((lb):(ub),:)) ;
                real(AllFramesFTrealign((lb):(ub),:)) ];
                %don't display rejects
                plotrealign((ub-lb+1):end,(MRS_struct.reject(:,ii).'==1))=min(plotrealign(:));
                imagesc(plotrealign);
                title('Cr Frequency, pre and post align');
                xlabel('time');
                 set(gca,'YTick',[1 CrFitRange CrFitRange+CrFitRange*(CrFitLimHigh-3.02)/(CrFitLimHigh-CrFitLimLow) CrFitRange*2]);
                 set(gca,'YTickLabel',[CrFitLimHigh CrFitLimLow 3.02 CrFitLimLow]);
                 %Add in labels for pre post
                 text(size(plotrealign,2)/18*17,0.4*size(plotrealign,1), 'PRE', 'Color',[1 1 1],'HorizontalAlignment','right');
                 text(size(plotrealign,2)/18*17,0.9*size(plotrealign,1), 'POST', 'Color',[1 1 1],'HorizontalAlignment','right');
             else
                 tmp = 'No realignment';
                 text(0,0.9, tmp, 'FontName', 'Courier');
             end

             %Bottom Right
             subplot(2,2,4);
             axis off;
             if strcmp(MRS_struct.vendor,'Siemens')
                 tmp = [ 'filename    : ' MRS_struct.gabafile{ii*2-1} ];
             else
                tmp = [ 'filename    : ' MRS_struct.gabafile{ii} ];
             end
             tmp = regexprep(tmp, '_','-');
             text(0,0.9, tmp, 'FontName', 'Helvetica','FontSize',12);
             tmp = [ 'Navg        : ' num2str(MRS_struct.Navg(ii)) ];
             text(0,0.8, tmp, 'FontName', 'Helvetica','FontSize',12);
             tmp = sprintf('Cr FWHM (Hz)   : %.2f', MRS_struct.CrFWHMHz(ii) );
             text(0,0.7, tmp, 'FontName', 'Helvetica','FontSize',12);
             %tmp = sprintf('FreqSTD (Hz): %.2f', MRS_struct.FreqStdevHz(ii));
             %text(0,0.6, tmp, 'FontName', 'Helvetica','FontSize',12);
             tmp = [ 'LB (Hz)     : ' num2str(MRS_struct.LB,2) ];
             text(0,0.5, tmp, 'FontName', 'Helvetica','FontSize',12);
             %tmp = [ 'Align/Reject: ' num2str(MRS_struct.FreqPhaseAlign) ];
             %text(0,0.5, tmp, 'FontName', 'Courier');
             tmp = [ 'Rejects     : '  num2str(sum(MRS_struct.reject(:,ii),1)) ];
             text(0,0.4, tmp, 'FontName', 'Courier');
             tmp = [ 'LoadVer     : ' MRS_struct.versionload ];
             text(0,0.3, tmp, 'FontName', 'Helvetica','FontSize',12);
    %         
              script_path=which('GannetLoad');
              % CJE update for GE
    %          Gannet_circle=[script_path(1:(end-12)) 'GANNET_circle.png'];
              Gannet_circle_white=[script_path(1:(end-13)) '/GANNET_circle_white.jpg'];
    %          A=imread(Gannet_circle);
              A2=imread(Gannet_circle_white);
              hax=axes('Position',[0.80, 0.05, 0.15, 0.15]);
              %set(gca,'Units','normalized');set(gca,'Position',[0.05 0.05 1.85 0.15]);
              image(A2);axis off; axis square;
              if strcmp(MRS_struct.vendor,'Siemens')
                  pfil_nopath = MRS_struct.gabafile{ii*2-1};
              else
                pfil_nopath = MRS_struct.gabafile{ii};
              end
              %for philips .data
              if(strcmpi(MRS_struct.vendor,'Philips_data'))
              fullpath = MRS_struct.gabafile{ii};
              %fullpath = regexprep(fullpath, '\./', '');   edit out NP   
              %fullpath = regexprep(fullpath, '/', '_'); edit out NP
              fullpath = regexprep(fullpath, '.data', '_data'); %NP see below
              fullpath = regexprep(fullpath, '\', '_');
              fullpath = regexprep(fullpath, '/', '_');
              %NP edit 02012013
              %Previous code somehow didn't run when running from hierarchical
              %folder (e.g. GABA_file = '.\name\MRI\raw.data) I got an error when Gannet tried to save the pdf for
              %.data file. E.g. ??? Error using ==> saveas at 115 Invalid or missing path: MRSload_output/.\7011-0124\MRI\raw_008.data.pdf
              %So it obviously didn't rewrite the path properly for the pdf here, but it IS important to get both folder and filename
              %as a lot of the .data files have similar names (e.g.
              %%raw_001.data). This change works for me for now, might not
              %%be most elegant
              
              end
              %  pfil_nopath = pfil_nopath( (length(pfil_nopath)-15) : (length(pfil_nopath)-9) );
              tmp = strfind(pfil_nopath,'/');
              tmp2 = strfind(pfil_nopath,'\');
              if(tmp)
                  lastslash=tmp(end);
              elseif (tmp2)
                  %maybe it's Windows...
                  lastslash=tmp2(end);
              else
                  % it's in the current dir...
                  lastslash=0;
              end
    %           
               if(strcmpi(MRS_struct.vendor,'Philips'))
                   tmp = strfind(pfil_nopath, '.sdat');
                   tmp1= strfind(pfil_nopath, '.SDAT');
                   if size(tmp,1)>size(tmp1,1)
                       dot7 = tmp(end); % just in case there's another .sdat somewhere else...
                   else
                       dot7 = tmp1(end); % just in case there's another .sdat somewhere else...
                   end
               elseif(strcmpi(MRS_struct.vendor,'GE'))
                  tmp = strfind(pfil_nopath, '.7');
                  dot7 = tmp(end); % just in case there's another .7 somewhere else...
              elseif(strcmpi(MRS_struct.vendor,'Philips_data'))
                  tmp = strfind(pfil_nopath, '.data');
                  dot7 = tmp(end); % just in case there's another .data somewhere else...
              elseif(strcmpi(MRS_struct.vendor,'Siemens'))
                  tmp = strfind(pfil_nopath, '.rda');
                  dot7 = tmp(end); % just in case there's another .rda somewhere else...
              elseif(strcmpi(MRS_struct.vendor,'Siemens_twix'))
                  tmp = strfind(pfil_nopath, '.dat');
                  dot7 = tmp(end); % just in case there's another .rda somewhere else...    
              end
              pfil_nopath = pfil_nopath( (lastslash+1) : (dot7-1) );
              %hax=axes('Position',[0.85, 0.05, 0.15, 0.15]);
              %set(gca,'Units','normalized');set(gca,'Position',[0.05 0.05 1.85 0.15]);
              %image(A2);axis off; axis square;
              % fix pdf output, where default is cm
              if sum(strcmp(listfonts,'Helvetica'))>0
               set(findall(h,'type','text'),'FontName','Helvetica');
               set(ha,'FontName','Helvetica');
               set(hb,'FontName','Helvetica');
               set(hc,'FontName','Helvetica');
              end
              set(gcf, 'PaperUnits', 'inches');
              set(gcf,'PaperSize',[11 8.5]);
              set(gcf,'PaperPosition',[0 0 11 8.5]);
              if(strcmpi(MRS_struct.vendor,'Philips_data'))
                  pdfname=[ 'MRSload_output/' fullpath '.pdf' ];
              else
                  pdfname=[ 'MRSload_output/' pfil_nopath  '.pdf' ];
              end
              saveas(h, pdfname);



                  %hax=axes('Position',[0.85, 0.05, 0.15, 0.15]);
         %          %set(gca,'Units','normalized');set(gca,'Position',[0.05 0.05 1.85 0.15]);

        %          
        %          
        %          
        %         
        %         
        %         
        %         
        %         
        %         











                %Save the structure?





          end%end of load-and-processing loop over datasets
end