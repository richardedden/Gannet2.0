function MRS_struct = GannetSegment(MRS_struct)

% Relies on SPM being installed
%
% Runs segmentation script if segmented images not present according to
% file convention of c1 c2 and c3 as prefixes on the anatomical image name
% for the GM, WM and CSF segmentations. If these files are present, they are
% loaded and used for the voxel segmentation

MRS_struct.version.segment = '170831';

if MRS_struct.p.PRIAM % deciding how many regions are there -- MGSaleh 2016
    vox = {MRS_struct.p.Vox};
else
    vox = {MRS_struct.p.Vox{1}};
end

% Set up SPM for batch processing
spm('defaults','fmri');
spm_jobman('initcfg');

for ii = 1:length(MRS_struct.metabfile)    
        
        % 1 - Take nifti from GannetCoRegister and segment it in SPM
        
        [T1dir, T1name, T1ext] = fileparts(MRS_struct.mask.T1image{ii});
        anatimage = MRS_struct.mask.T1image{ii};
        
        % Check to see if segmentation already done - if not, do it
        % Check which SPM version is installed and segment accordingly
        tmp = [T1dir '/c1' T1name T1ext];
        if ~exist(tmp,'file')
            spmversion = fileparts(which('spm'));
            if strcmpi(spmversion(end-4:end),'spm12')
                CallSPM12segmentation(anatimage);
            else
                CallSPM8segmentation(anatimage);
            end
        end
        
        % 2 - Determine GM, WM and CSF fractions for each voxel
        
        if strcmp(T1dir,'')
            T1dir='.';
        end
        
        GM  = [T1dir '/c1' T1name T1ext];
        WM  = [T1dir '/c2' T1name T1ext];
        CSF = [T1dir '/c3' T1name T1ext];
        
        GMvol  = spm_vol(GM);
        WMvol  = spm_vol(WM);
        CSFvol = spm_vol(CSF);        
        
        % Loop over voxels if PRIAM
        for kk = 1:length(vox)
            
            voxmaskvol = spm_vol(cell2mat(MRS_struct.mask.outfile(ii)));
            
            % GM
            O_GMvox.fname = [T1dir '/c1' T1name '_GM_mask.nii'];
            O_GMvox.descrip = 'GMmasked_MRS_Voxel_Mask';
            O_GMvox.dim = voxmaskvol.dim;
            O_GMvox.dt = voxmaskvol.dt;
            O_GMvox.mat = voxmaskvol.mat;
            GM_voxmask_vol = GMvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
            O_GMvox = spm_write_vol(O_GMvox, GM_voxmask_vol);
            
            % WM
            O_WMvox.fname = [T1dir '/c2' T1name '_WM_mask.nii'];
            O_WMvox.descrip = 'WMmasked_MRS_Voxel_Mask';
            O_WMvox.dim = voxmaskvol.dim;
            O_WMvox.dt = voxmaskvol.dt;
            O_WMvox.mat = voxmaskvol.mat;
            WM_voxmask_vol = WMvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
            O_WMvox = spm_write_vol(O_WMvox, WM_voxmask_vol);
            
            % CSF
            O_CSFvox.fname = [T1dir '/c3' T1name '_CSF_mask.nii'];
            O_CSFvox.descrip = 'CSFmasked_MRS_Voxel_Mask';
            O_CSFvox.dim = voxmaskvol.dim;
            O_CSFvox.dt = voxmaskvol.dt;
            O_CSFvox.mat = voxmaskvol.mat;
            CSF_voxmask_vol = CSFvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
            O_CSFvox = spm_write_vol(O_CSFvox, CSF_voxmask_vol);
            
            % ***NB***
            % For a subject with multiple voxels, the segmented voxels will
            % get overwritten
            
            % 3 - Calculate an adjusted gabaiu and output it to the structure
            
            gmsum = sum(sum(sum(O_GMvox.private.dat(:,:,:))));
            wmsum = sum(sum(sum(O_WMvox.private.dat(:,:,:))));
            csfsum = sum(sum(sum(O_CSFvox.private.dat(:,:,:))));
            
            gmfra = gmsum/(gmsum+wmsum+csfsum);
            wmfra = wmsum/(gmsum+wmsum+csfsum);
            csffra = csfsum/(gmsum+wmsum+csfsum);
            
            tissuefra = gmfra+wmfra;
            
            MRS_struct.out.(vox{kk}).tissue.GMfra(ii) = gmfra;
            MRS_struct.out.(vox{kk}).tissue.WMfra(ii) = wmfra;
            MRS_struct.out.(vox{kk}).tissue.CSFfra(ii) = csffra;
            
            % Correction of institutional units only feasible if water-scaling is
            % performed, skip otherwise (GO 07/13/2017)
            % MM (170831): Loop over voxels and metabolites
            if isfield(MRS_struct,'waterfile')                
                if MRS_struct.p.HERMES
                    target = {MRS_struct.p.target, MRS_struct.p.target2};
                else
                    target = {MRS_struct.p.target};
                end
                
                for trg = 1:length(target)
                    if strcmp(target{trg},'GABAGlx')
                        MRS_struct.out.(vox{kk}).GABA.ConcIU_TissCorr(ii) = ...
                            MRS_struct.out.(vox{kk}).GABA.ConcIU(ii) / tissuefra;
                        MRS_struct.out.(vox{kk}).Glx.ConcIU_TissCorr(ii) = ...
                            MRS_struct.out.(vox{kk}).Glx.ConcIU(ii) / tissuefra;
                    else
                        MRS_struct.out.(vox{kk}).(target{trg}).ConcIU_TissCorr(ii) = ...
                            MRS_struct.out.(vox{kk}).(target{trg}).ConcIU(ii) / tissuefra;
                    end
                end
            end
            
            % 4 - Build output
            
            if ishandle(104)
                clf(104) % MM (170831)
            end
            h = figure(104);
            % MM (170831): Open figure in center of screen
            scr_sz = get(0, 'ScreenSize');
            fig_w = 1000;
            fig_h = 707;
            set(h, 'Position', [(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
            set(h,'Color',[1 1 1]);
            figTitle = 'GannetSegment Output';
            set(gcf,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
            
            % Voxel co-registration
            subplot(2,2,1); % a subplot of the voxel on the brain
            imagesc(squeeze(MRS_struct.mask.img(ii,:,1:round(size(MRS_struct.mask.img,3)/3))));
            colormap('gray');
            caxis([0 0.5])
            axis equal;
            axis tight;
            axis off;
            
            subplot(2,2,3); % replot of GABA fit spec
            z=abs(MRS_struct.spec.freq-4.1);
            lowerbound=find(min(z)==z);
            z=abs(MRS_struct.spec.freq-2.79);
            upperbound=find(min(z)==z);
            freqbounds=lowerbound:upperbound;
            freq=MRS_struct.spec.freq(freqbounds);
            plot(MRS_struct.spec.freq, real(MRS_struct.spec.(vox{kk}).GABAGlx.diff(ii,:)), ...
                'k', freq, GABAGlxModel(MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:),freq), 'r');
            
            zz=abs(MRS_struct.spec.freq-3.6);
            Glx_right=find(min(zz)==zz);
            %zz=abs(MRS_struct.spec.freq-3.3);
            %GABA_left=find(min(zz)==zz);
            zz=abs(MRS_struct.spec.freq-2.8);
            GABA_right=find(min(zz)==zz);
            %specbaseline = mean(real(SpectraToPlot(1,GABA_right:GABA_left)),2);
            gabaheight = max(abs(real(MRS_struct.spec.(vox{kk}).GABAGlx.diff(ii,Glx_right:GABA_right))),[],2);
            gabaheight = mean(gabaheight);
            
            yaxismax = 2*gabaheight; % top spec + 2* height of gaba
            yaxismin = -2*gabaheight; % extend 2* gaba heights below zero
            if yaxismax < yaxismin
                dummy=yaxismin;
                yaxismin=yaxismax;
                yaxismax=dummy;
            end
            axis([0 5 yaxismin yaxismax])
            set(gca,'YTick',[]);
            set(gca,'XLim',[0 4.5]);
            set(gca,'XDir','reverse');
            
            subplot(2,2,2)  % output results
            axis off;
            
            % Print correction of institutional units only feasible if water-scaling is
            % performed, skip otherwise (GO 07/13/2017)
            if isfield(MRS_struct,'waterfile')
                tmp = ['GABAconc(iu) tissue corr:  ' num2str(MRS_struct.out.(vox{kk}).GABA.ConcIU_TissCorr(ii))];
                text(0, 0.87, tmp, 'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'top',...
                    'FontName', 'Helvetica','FontSize',13);
            end
            
            tmp = ['GM voxel fraction:  ' num2str(MRS_struct.out.(vox{kk}).tissue.GMfra(ii))];
            text(0,0.75, tmp, 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top',...
                'FontName', 'Helvetica','FontSize',13);
            tmp = ['WM voxel fraction:  ' num2str(MRS_struct.out.(vox{kk}).tissue.WMfra(ii))];
            text(0,0.63, tmp, 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top',...
                'FontName', 'Helvetica','FontSize',13);
            tmp = ['CSF voxel fraction:  ' num2str(MRS_struct.out.(vox{kk}).tissue.CSFfra(ii))];
            text(0,0.51, tmp, 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top',...
                'FontName', 'Helvetica','FontSize',13);
            
            C = MRS_struct.metabfile{ii};
            if size(C,2) > 30
                [~,y] = fileparts(C);
            else
                y = C;
            end
            tmp = ['MRS data:  ' y];
            tmp = regexprep(tmp, '_','-');
            text(0,0.39, tmp, 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top',...
                'FontName', 'Helvetica','FontSize',13);
            
            D = MRS_struct.mask.T1image{ii} ;
            if size(D,2) > 30
                [~,y] = fileparts(D);
            else
                y = D;
            end
            tmp = ['Anatomical image:  ' y];
            tmp = regexprep(tmp, '_','-');
            text(0,0.27, tmp, 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top',...
                'FontName', 'Helvetica','FontSize',13);
            
            tmp = ['SegmentVer:  ' MRS_struct.version.segment];
            text(0,0.15, tmp, 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top',...
                'FontName', 'Helvetica','FontSize',13);
            
            subplot(2,2,4);
            axis off;
            script_path=which('GannetFit');
            Gannet_logo=[script_path(1:(end-12)) '/Gannet3_logo.png'];
            A2=imread(Gannet_logo,'png','BackgroundColor',[1 1 1]);
            axes('Position',[0.80, 0.05, 0.15, 0.15]);
            image(A2); axis off; axis square;
            
            %%%% Save PDF %%%%%
            pdfdirname = './GannetSegment_output'; % MM (170831)
            pfil_nopath = MRS_struct.metabfile{ii};
            
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
            if strcmpi(MRS_struct.p.vendor,'Philips')
                tmp = strfind(pfil_nopath, '.sdat');
                tmp1= strfind(pfil_nopath, '.SDAT');
                if size(tmp,1)>size(tmp1,1)
                    dot7 = tmp(end); % just in case there's another .sdat somewhere else...
                else
                    dot7 = tmp1(end); % just in case there's another .sdat somewhere else...
                end
            elseif(strcmpi(MRS_struct.p.vendor,'GE'))
                tmp = strfind(pfil_nopath, '.7');
                dot7 = tmp(end); % just in case there's another .7 somewhere else...
            elseif(strcmpi(MRS_struct.p.vendor,'Philips_data'))
                tmp = strfind(pfil_nopath, '.data');
                dot7 = tmp(end); % just in case there's another .data somewhere else...
            elseif(strcmpi(MRS_struct.p.vendor,'Siemens'))
                tmp = strfind(pfil_nopath, '.rda');
                dot7 = tmp(end); % just in case there's another .rda somewhere else...
            elseif(strcmpi(MRS_struct.p.vendor,'Siemens_twix'))
                tmp = strfind(pfil_nopath, '.dat');
                dot7 = tmp(end); % just in case there's another .dat somewhere else...
            end
            pfil_nopath = pfil_nopath(lastslash+1:dot7-1);
            
            %Save PDF output
            set(gcf,'PaperUnits','inches');
            set(gcf,'PaperSize',[11 8.5]);
            set(gcf,'PaperPosition',[0 0 11 8.5]);
            if(strcmpi(MRS_struct.p.vendor,'Philips_data'))
                pdfname = [pdfdirname '/' fullpath '_segment.pdf'];
            else
                pdfname = [pdfdirname '/' pfil_nopath  '_segment.pdf'];
            end
            if ~exist(pdfdirname,'dir')
                mkdir(pdfdirname)
            end
            saveas(gcf, pdfname);
            
        end
        
end

end

