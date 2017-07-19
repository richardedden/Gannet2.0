function [MRS_struct] =  GannetSegment(MRS_struct)

% relies on spm being installed
%
% Runs segmentation script if segmented images not present according to
% file convention of c1 c2 and c3 as prefixes on the anatomical image name
% for the GM WM and CSF segmentations. If these files are present, they are
% loaded and used for the voxel segmentation

MRS_struct.out.tissue.version= '20170713';

for ii = 1:MRS_struct.ii
    
[T1dir, T1name, T1ext] = fileparts(MRS_struct.mask.T1image{ii});

%1 - take nifti from GannetCoRegister and segment it in spm

anatimage = MRS_struct.mask.T1image{ii};

%check to see if segmentation done - if its not done, do it

tmp = [T1dir '/c1' T1name T1ext];
if exist(tmp,'var') ~= 2
    CallSPMsegmentation(anatimage);
    garbage = [T1dir '/c4' T1name T1ext];
    delete(garbage)
    garbage = [T1dir '/c5' T1name T1ext];
    delete(garbage)
end

%2 - determing GM,WM and CSF fractions for each voxel

if strcmp(T1dir,'')
    T1dir='.';
end

GM = [T1dir '/c1' T1name T1ext];
WM = [T1dir '/c2' T1name T1ext];
CSF = [T1dir '/c3' T1name T1ext];

GMvol = spm_vol(GM);
WMvol = spm_vol(WM);
CSFvol = spm_vol(CSF);

voxmaskvol = spm_vol(cell2mat(MRS_struct.mask.outfile(ii)));

%GM
O_GMvox.fname = [T1dir '/c1' T1name '_GM.nii'];
O_GMvox.descrip='GMmasked_MRS_Voxel_Mask';
O_GMvox.dim = voxmaskvol.dim;
O_GMvox.dt = voxmaskvol.dt; 
O_GMvox.mat = voxmaskvol.mat;
GM_voxmask_vol = GMvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
O_GMvox = spm_write_vol(O_GMvox, GM_voxmask_vol);

%WM
O_WMvox.fname = [T1dir '/c2' T1name '_WM.nii'];
O_WMvox.descrip='WMmasked_MRS_Voxel_Mask';
O_WMvox.dim = voxmaskvol.dim;
O_WMvox.dt = voxmaskvol.dt; 
O_WMvox.mat = voxmaskvol.mat;
WM_voxmask_vol = WMvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
O_WMvox = spm_write_vol(O_WMvox, WM_voxmask_vol);

%CSF
O_CSFvox.fname = [T1dir '/c3' T1name '_CSF.nii'];
O_CSFvox.descrip='CSFmasked_MRS_Voxel_Mask';
O_CSFvox.dim = voxmaskvol.dim;
O_CSFvox.dt = voxmaskvol.dt; 
O_CSFvox.mat = voxmaskvol.mat;
CSF_voxmask_vol = CSFvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
O_CSFvox = spm_write_vol(O_CSFvox, CSF_voxmask_vol);

% ***NB*** for a subject with multiple voxel, the segmented voxels will get
%overwritten


%3 - calculate an adjusted gabaiu and output it to the structure

gmsum = sum(sum(sum(O_GMvox.private.dat(:,:,:))));
wmsum = sum(sum(sum(O_WMvox.private.dat(:,:,:))));
csfsum = sum(sum(sum(O_CSFvox.private.dat(:,:,:))));

gmfra = gmsum/(gmsum+wmsum+csfsum);
wmfra = wmsum/(gmsum+wmsum+csfsum);
csffra = csfsum/(gmsum+wmsum+csfsum);

tissuefra = gmfra+wmfra;

MRS_struct.out.tissue.GMfra(ii) = gmfra;
MRS_struct.out.tissue.WMfra(ii) = wmfra;
MRS_struct.out.tissue.CSFfra(ii) = csffra;

% Correction of institutional units only feasible if water-scaling is 
% performed, skip otherwise (GO 07/13/2017)
if isfield(MRS_struct,'waterfile')
    MRS_struct.out.GABAconciuTissCorr(ii) = MRS_struct.out.GABAconciu(ii)./tissuefra;
end

%4 - build output


fignum = 104;
if(ishandle(fignum))
close(fignum)
end
h=figure(fignum);

set(h, 'Position', [100, 100, 1000, 707]);
set(h,'Color',[1 1 1]);
figTitle = 'GannetSegment Output';
set(gcf,'Name',figTitle,'Tag',figTitle, 'NumberTitle','off');
% GABA plot

subplot(2, 2, 1);    % a subplot of the voxel on the brain

% input=MRS_struct.mask.img(ii,:,1:round(size(MRS_struct.mask.img,3)/3));
% size(input)
imagesc(squeeze(MRS_struct.mask.img(ii,:,1:round(size(MRS_struct.mask.img,3)/3))));
colormap('gray');
caxis([0 0.5])
axis equal;
axis tight;
axis off;
 

subplot(2, 2, 3);    % replot of GABA fit spec
z=abs(MRS_struct.spec.freq-3.55);
lowerbound=find(min(z)==z);
z=abs(MRS_struct.spec.freq-2.79);
upperbound=find(min(z)==z);
freqbounds=lowerbound:upperbound;
freq=MRS_struct.spec.freq(1,freqbounds);
plot( ...
    real(MRS_struct.spec.freq(1,:)),real(MRS_struct.spec.diff(ii,:)), ...
    'k',freq,GaussModel(MRS_struct.out.GABAModelFit(ii,:),freq),'r');  % this part may be broken
        
zz=abs(MRS_struct.spec.freq-3.6);
Glx_right=find(min(zz)==zz);
%zz=abs(MRS_struct.spec.freq-3.3);
%GABA_left=find(min(zz)==zz);
zz=abs(MRS_struct.spec.freq-2.8);
GABA_right=find(min(zz)==zz);
%specbaseline = mean(real(SpectraToPlot(1,GABA_right:GABA_left)),2);
gabaheight = abs(max(MRS_struct.spec.diff(1,Glx_right:GABA_right),[],2));
gabaheight = mean(gabaheight);

yaxismax =  2 *gabaheight; % top spec + 2* height of gaba
yaxismin =  -2* gabaheight; % extend 2* gaba heights below zero
if (yaxismax<yaxismin)
    dummy=yaxismin;
    yaxismin=yaxismax;
    yaxismax=dummy;
end
axis([0 5  yaxismin yaxismax])
set(gca,'YTick',[]);
set(gca,'XLim',[0 4.5]);
set(gca,'XDir','reverse');
%    set(gca, 'Box', 'off');

subplot(2,2,2)  % output results
axis off;

% Print correction of institutional units only feasible if water-scaling is 
% performed, skip otherwise (GO 07/13/2017)
if isfield(MRS_struct,'waterfile')
    tmp = ['GABAconc(iu) tissue corr:  ' num2str(MRS_struct.out.GABAconciuTissCorr(ii))];
    text(0, 0.87, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
end

tmp = ['Voxel fraction GM:  ' num2str(MRS_struct.out.tissue.GMfra(ii))];
text(0,0.75, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
tmp = ['Voxel fraction WM:  ' num2str(MRS_struct.out.tissue.WMfra(ii))];
text(0,0.63, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
tmp = ['Voxel fraction CSF:  ' num2str(MRS_struct.out.tissue.CSFfra(ii))];
text(0,0.51, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);

C = MRS_struct.gabafile{ii}; 
if size(C,2) >30
    [x, y, z] = fileparts(C);
else
    y = MRS_struct.gabafile(ii);
end
tmp = ['Spectra:  ' y];
tmp = regexprep(tmp, '_','-');
text(0,0.39, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        
D = MRS_struct.mask.T1image{ii} ;
if size(D,2) >30
    [x, y, z] = fileparts(D);
else
    y = MRS_struct.gabafile(ii);
end
tmp = ['Anatomical image:  ' MRS_struct.mask.T1image(ii,:)];
tmp = regexprep(tmp, '_','-');
text(0,0.27, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        
        
        
tmp = ['Ver(Segment):  ' (MRS_struct.out.tissue.version)];
text(0,0.15, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);

            
subplot(2,2,4)
    axis off
    script_path=which('GannetFit');
    Gannet_circle_white=[script_path(1:(end-12)) '/GANNET_circle_white.jpg'];
    A_2=imread(Gannet_circle_white);
    axes('Position',[0.80, 0.05, 0.15, 0.15]);
    image(A_2);axis off; axis square;


    
               %%%%  Save EPS %%%%%
           

epsdirname = [ './MRSSegment_' datestr(clock,'yymmdd') ];

pfil_nopath = MRS_struct.gabafile{ii};

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
    if(strcmpi(MRS_struct.p.vendor,'Philips'))
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
    pfil_nopath = pfil_nopath( (lastslash+1) : (dot7-1) );

    
    
    %Save pdf output
    
set(gcf, 'PaperUnits', 'inches');
    set(gcf,'PaperSize',[11 8.5]);
    set(gcf,'PaperPosition',[0 0 11 8.5]);
    if(strcmpi(MRS_struct.p.vendor,'Philips_data'))
        pdfname=[ epsdirname '/' fullpath '.pdf' ];
    else
        pdfname=[ epsdirname '/' pfil_nopath  '.pdf' ];
    end
    %epsdirname
    if(exist(epsdirname,'dir') ~= 7)
        %epsdirname
        mkdir(epsdirname)
    end
    saveas(gcf, pdfname);
%     
% if(ii==MRS_struct.ii)
%     if(isfield(MRS_struct.p, 'mat') == 1)
%        matname =[ epsdirname '/' 'MRS_struct' '.mat' ];
%        save(matname,'MRS_struct'); 
%     end
% end
% 
%         
  

end

end

