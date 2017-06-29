function MRS_struct = GannetQuantify(MRS_struct)

% objective: include segmentation of GM, WM and CSF components to integrate
% T1 and T2 of water in these tissues properly. 

% need to bring in MRS_struct 
% need to re-develop the GABA_iu quantify equation - see notebook
% need to make sure use proper T1 and T2 for water in GM, WM and CSF - see
% notebook for values

% then integrate this into MRS_struct

TR = MRS_struct.p.TR/1000;
TE = MRS_struct.p.TE/1000;
meanGMfra = mean(MRS_struct.out.tissue.GMfra);  % fraction of voxel that is GM for voxel fraction normalization - move to preinitialize
meanWMfra = mean(MRS_struct.out.tissue.WMfra);
cWM = 1;  % concentration of GABA in WM
cGM = 2;  % concentration of GABA in GM
          % note: cWM/cGM is term used ratio not absolute value is
          % important piece
alpha = cWM/cGM;
% Constants
% from Wansapura 1999  JMRI; 9:531 (1999)
%        T1          T2
% WM   832 +/- 10  79.2 +/- 0.6
% GM  1331 +/- 13  110 +/- 2
%
% from Lu, JMRI; 2005; 22: 13
% CSF T1 = 3817 +/- 424msec - but state may underestimated and that 4300ms 
% is likely more accurate - but the reference is to an ISMRM 2001 abstract
% MacKay (last author) 2006 ISMRM abstract has T1 CSF = 3300 ms 
% CSF T2 = 503.0 +/- 64.3 Piechnik MRM 2009; 61: 579
% However, other values from Stanisz et al:
% CPMG for T2, IR for T1
% T2GM = 99 +/ 7, lit: 71+/- 10 (27)
% T1GM = 1820 +/- 114, lit 1470 +/- 50 (29)
% T2WM = 69 +/-3 lit 56 +/- 4 (27)
% T1WM = 1084 +/- 45 lit 1110 +/- 45 (29)

% Constants (s)
T1w_WM = 0.832;
T2w_WM = 0.0792;
T1w_GM = 1.331;
T2w_GM = 0.110;
T1w_CSF = 3.817;
T2w_CSF = 0.503;
% determine concentration of water in GM WM and CSF
% Gasparovic et al, MRM 2006; 55:1219 uses relative densities, ref to Ernst
% fGM = 0.78
% fWM = 0.65 
% fCSH = 0.97
% such that 
% concw_GM = 0.78 * 55.51 mol/kg = 43.30
% concw_WM = 0.65 * 55.51 mol/kg = 36.08
% concw_CSF = 0.97 * 55.51 mol/kg = 53.84



T1_GABA = 0.80 ; % "empirically determined"...! Gives same values as RE's spreadsheet
% ... and consistent with Cr-CH2 T1 of 0.8 (Traber, 2004)
%Not yet putting in measured GABA T1, but it is in the pipeline - 1.35ish

%T2_GABA = 0.13; % from occipital Cr-CH2, Traber 2004
T2_GABA = 0.088; % from JMRI paper 2011 Eden et al.

concw_GM = 43.30*1000;
concw_WM = 36.08*1000;
concw_CSF = 53.84*1000;


for ii = 1: MRS_struct.ii


EditingEff = 0.5; % GannetFit, Gannet1
MM = 0.45;        % GannetFit, Gannet1
fracGM = MRS_struct.out.tissue.GMfra(ii);
fracWM = MRS_struct.out.tissue.WMfra(ii);
fracCSF = MRS_struct.out.tissue.CSFfra(ii);

CorrFactor = (meanGMfra + alpha* meanWMfra) /( (fracGM + alpha*fracWM)*(meanGMfra + meanWMfra)) ;

MRS_struct.Quantify.QuantGABA_iu(ii) = (MRS_struct.out.GABAArea(ii)  ./  MRS_struct.out.WaterArea(ii))  ...
    * MM / EditingEff *  ...
    (  ...
    fracGM * concw_GM * (1-exp(-TR/T1w_GM)) * (exp(-TE/T2w_GM))/( (1-exp(-TR/T1_GABA)) * (exp(-TE/T2_GABA)))  ...
    +  ...
    fracWM *concw_WM * (1-exp(-TR/T1w_WM)) * (exp(-TE/T2w_WM))/( (1-exp(-TR/T1_GABA)) * (exp(-TE/T2_GABA)))  ...
    +      ...
    fracCSF *concw_CSF * (1-exp(-TR/T1w_CSF)) * (exp(-TE/T2w_CSF))/( (1-exp(-TR/T1_GABA)) * (exp(-TE/T2_GABA)))  ...
    ) ;

MRS_struct.Quantify.QuantCorrGABA_iu(ii) = MRS_struct.Quantify.QuantGABA_iu(ii) / (fracGM +alpha*fracWM);

MRS_struct.Quantify.QuantNormTissCorrGABA_iu(ii) =  MRS_struct.Quantify.QuantGABA_iu(ii) *  CorrFactor;
    
% MRS_struct.Quantify.tissalpha = MRS_struct.out.GABAconciu * CorrFactor;
% probably make a new output - combination of GannetSeg and what did
% here...


% save .mat 

epsdirname = [ './MRSQuantify_' datestr(clock,'yymmdd') ];
    if(exist(epsdirname,'dir') ~= 7)
        %epsdirname
        mkdir(epsdirname)
    end

if(isfield(MRS_struct.p, 'mat') == 1)
       matname =[ epsdirname '/' 'MRS_struct' '.mat' ];
       save(matname,'MRS_struct'); 
end




% build output pdf summary

fignum = 501;
if(ishandle(fignum))
close(fignum)
end
h=figure(fignum);

set(h, 'Position', [100, 100, 1000, 707]);
set(h,'Color',[1 1 1]);
figTitle = 'GannetQuantify Output';
set(gcf,'Name',figTitle,'Tag',figTitle, 'NumberTitle','off');
% GABA plot

% a subplot of the voxel on the brain
% replot of GABA fit spec

% OUTPUTS:
% data files - .data and .mask
% GABA CSF corr
% GABA reln cnsts
% GABA alpha corr
% GABA normalize voxel alpha corr

 
h=subplot(2, 2, 4);    

% input=MRS_struct.mask.img(ii,:,1:round(size(MRS_struct.mask.img,3)/3));
% size(input)
imagesc(squeeze(MRS_struct.mask.img(ii,:,1:round(size(MRS_struct.mask.img,3)/3))));
colormap('gray');
caxis([0 0.5])
axis equal;
axis tight;
axis off;
 

h=subplot(2, 2, 1);   
    z=abs(MRS_struct.spec.freq-3.55);
    lowerbound=find(min(z)==z);
    z=abs(MRS_struct.spec.freq-2.79);
    upperbound=find(min(z)==z);
    freqbounds=lowerbound:upperbound;
    freq=MRS_struct.spec.freq(1,freqbounds);
    plot( ...
        real(MRS_struct.spec.freq(1,:)),real(MRS_struct.spec.diff(ii,:)), ...
        'k',freq,GaussModel_area(MRS_struct.out.GABAModelFit(ii,:),freq),'r');  % this part may be broken
        
zz=abs(MRS_struct.spec.freq-3.6);
Glx_right=find(min(zz)==zz);
zz=abs(MRS_struct.spec.freq-3.3);
GABA_left=find(min(zz)==zz);
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

% tmp = ['GABAconc(iu) tissue corr (CSF corrected):  ' num2str(MRS_struct.out.GABAconciuTissCorr(ii))];
%     text(0, 0.87, tmp, 'HorizontalAlignment', 'left', ...
%             'VerticalAlignment', 'top',...
%             'FontName', 'Helvetica','FontSize',13);

tmp = ['GABAconc(iu) with tissue relaxation and vis:  ' num2str(MRS_struct.Quantify.QuantGABA_iu(ii))];
    text(0, 0.77, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);

tmp = ['GABAconc(iu) alpha corrected:  ' num2str(MRS_struct.Quantify.QuantCorrGABA_iu(ii))];
    text(0, 0.67, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
tmp = ['    (uses tissue specific relaxation and visibility constants)' ];
    text(0, 0.6, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
tmp = ['GABAconc(iu) a-corrected, ave voxel norm:  ' num2str(MRS_struct.Quantify.QuantNormTissCorrGABA_iu(ii))];
    text(0, 0.5, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
tmp = ['    (uses tissue specific relaxation and visibility constants)' ];
    text(0, 0.43, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);


        
MRS_struct.Quantify.QuantCorrGABA_iu(ii) = MRS_struct.Quantify.QuantGABA_iu(ii) / (fracGM +alpha*fracWM);

MRS_struct.Quantify.QuantNormTissCorrGABA_iu(ii) =  MRS_struct.Quantify.QuantGABA_iu(ii) *  CorrFactor ;


C = MRS_struct.gabafile{ii}; 
if size(C,2) >30
    [x y z] = fileparts(C);
else
    y = MRS_struct.gabafile(ii);
end
tmp = ['Spectra:  ' y];
tmp = regexprep(tmp, '_','-');
text(0,0.25, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        
D = MRS_struct.mask.T1image{ii} ;
if size(D,2) >30
    [x y z] = fileparts(D);
else
    y = MRS_struct.gabafile(ii);
end
tmp = ['Anat image:  ' MRS_struct.mask.T1image(ii,:)];
tmp = regexprep(tmp, '_','-');
text(0,0.15, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);
        
        
        
tmp = ['Ver(Quantify):  ADH1' ];
text(0,0.0, tmp, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top',...
            'FontName', 'Helvetica','FontSize',13);

subplot(2,2,3)        
    rejectframesplot = (1./MRS_struct.out.reject(:,ii).') .*  MRS_struct.fids.waterfreq(ii,:);
    plot([1:size(MRS_struct.fids.data,2)], ...
        MRS_struct.fids.waterfreq(ii,:)', '-', ...
        [1:size(MRS_struct.fids.data,2)], rejectframesplot, 'ro');
    set(gca,'XLim',[0 size(MRS_struct.fids.data,2)]);
    xlabel('time'); ylabel('\omega_0');
    title('Water Frequency, ppm');

        
            
subplot(2,2,4)
    axis off
    script_path=which('GannetFit');
    Gannet_circle_white=[script_path(1:(end-12)) '/GANNET_circle_white.jpg'];
    A_2=imread(Gannet_circle_white);
    hax=axes('Position',[0.80, 0.05, 0.15, 0.15]);
    image(A_2);axis off; axis square;


    
               %%%%  Save EPS %%%%%
           
pfil_nopath = MRS_struct.gabafile{ii};
    
fullpath = MRS_struct.gabafile{ii};

    
if(strcmp(fullpath(1:2) , '..'))
fullpath = fullpath(4:end);
end

if(strcmpi(MRS_struct.p.vendor,'Philips_data'))
     fullpath = MRS_struct.gabafile{ii};    
    if(strcmp(fullpath(1:2) , '..'))
            fullpath = fullpath(4:end);
        end

   
    %         fullpath = regexprep(fullpath, '\./', ''); NP edit out.
    %         see below
    %         fullpath = regexprep(fullpath, '/', '_');
    fullpath = regexprep(fullpath, '.data', '_data');
    fullpath = regexprep(fullpath, '\', '_');
    fullpath = regexprep(fullpath, '/', '_');
    %NP edit 02012013
    %Previous code somehow didn't run when running from hierarchical
    %folder (e.g. GABA_file = '.\name\MRI\raw.data) I got an error when Gannet tried to save the pdf for
    %.data file. E.g. ??? Error using ==> saveas at 115 Invalid or missing path: ./MRSfit_140102/.\7011-0124\MRI\raw_008.data.pdf
    %So it obviously didn't rewrite the path properly for the pdf here, but it IS important to get both folder and filename
    %as a lot of the .data files have similar names (e.g.
    %%raw_001.data). This change works for me for now, might not
    %%be most elegant

end


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
    elseif(strcmpi(MRS_struct.p.vendor,'Philips_data'))  % make this be sdat
        tmp = strfind(pfil_nopath, '.data');
        dot7 = tmp(end); % just in case there's another .data somewhere else...
    end
    pfil_nopath = pfil_nopath( (lastslash+1) : (dot7-1) );

    
    
    %Save pdf output
    
    %Save pdf output
%
    
    
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
        

end
