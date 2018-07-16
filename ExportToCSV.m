function ExportToCSV(MRS_struct, target, kk, module)

round2 = @(x) round(x*1e3)/1e3;
vox = MRS_struct.p.Vox{kk};

for trg = 1:length(target)
    
    if ~strcmpi(target{trg},'GABAGlx')
        metab = target{trg};
    end
    
    if strcmp(MRS_struct.p.vendor,'Siemens_rda')
        out.matlabVersion = repmat(version('-release'), [length(MRS_struct.metabfile)/2 1]);
        out.GannetVersion = repmat('3', [length(MRS_struct.metabfile)/2 1]);
        out.date          = repmat(date, [length(MRS_struct.metabfile)/2 1]);
    else
        out.matlabVersion = repmat(version('-release'), [length(MRS_struct.metabfile) 1]);
        out.GannetVersion = repmat('3', [length(MRS_struct.metabfile) 1]);
        out.date          = repmat(date, [length(MRS_struct.metabfile) 1]);
    end
    
    %%% 1. Extract data from GannetFit %%%
    
    if strcmp(MRS_struct.p.vendor,'Siemens_rda')
        out.filename = MRS_struct.metabfile(1:2:end);
    else
        out.filename = MRS_struct.metabfile(:);
    end
    out.AvgDeltaF0 = MRS_struct.out.AvgDeltaF0(:);
    
    if strcmpi(target{trg},'GABAGlx')
        if strcmp(MRS_struct.p.Reference_compound,'H2O')
            out.GABA.ConcIU      = MRS_struct.out.(vox).GABA.ConcIU(:);
            out.GABA.ConcCr      = MRS_struct.out.(vox).GABA.ConcCr(:);
            out.GABA.FitError_W  = MRS_struct.out.(vox).GABA.FitError_W(:);
            out.GABA.FitError_Cr = MRS_struct.out.(vox).GABA.FitError_Cr(:);
            out.GABA.FWHM        = MRS_struct.out.(vox).GABA.FWHM(:);
            out.GABA.SNR         = MRS_struct.out.(vox).GABA.SNR(:);
            out.Glx.ConcIU       = MRS_struct.out.(vox).Glx.ConcIU(:);
            out.Glx.ConcCr       = MRS_struct.out.(vox).Glx.ConcCr(:);
            out.Glx.FitError_W   = MRS_struct.out.(vox).Glx.FitError_W(:);
            out.Glx.FitError_Cr  = MRS_struct.out.(vox).Glx.FitError_Cr(:);
            out.Glx.FWHM         = MRS_struct.out.(vox).Glx.FWHM(:);
            out.Glx.SNR          = MRS_struct.out.(vox).Glx.SNR(:);
        else
            out.GABA.ConcCr      = MRS_struct.out.(vox).GABA.ConcCr(:);
            out.GABA.FitError_Cr = MRS_struct.out.(vox).GABA.FitError_Cr(:);
            out.GABA.FWHM        = MRS_struct.out.(vox).GABA.FWHM(:);
            out.GABA.SNR         = MRS_struct.out.(vox).GABA.SNR(:);
            out.Glx.ConcCr       = MRS_struct.out.(vox).Glx.ConcCr(:);
            out.Glx.FitError_Cr  = MRS_struct.out.(vox).Glx.FitError_Cr(:);
            out.Glx.FWHM         = MRS_struct.out.(vox).Glx.FWHM(:);
            out.Glx.SNR          = MRS_struct.out.(vox).Glx.SNR(:);
        end
    else
        if strcmp(MRS_struct.p.Reference_compound,'H2O')
            out.(metab).ConcIU      = MRS_struct.out.(vox).(metab).ConcIU(:);
            out.(metab).ConcCr      = MRS_struct.out.(vox).(metab).ConcCr(:);
            out.(metab).FitError_W  = MRS_struct.out.(vox).(metab).FitError_W(:);
            out.(metab).FitError_Cr = MRS_struct.out.(vox).(metab).FitError_Cr(:);
            out.(metab).FWHM        = MRS_struct.out.(vox).(metab).FWHM(:);
            out.(metab).SNR         = MRS_struct.out.(vox).(metab).SNR(:);
        else
            out.(metab).ConcCr      = MRS_struct.out.(vox).(metab).ConcCr(:);
            out.(metab).FitError_Cr = MRS_struct.out.(vox).(metab).FitError_Cr(:);
            out.(metab).FWHM        = MRS_struct.out.(vox).(metab).FWHM(:);
            out.(metab).SNR         = MRS_struct.out.(vox).(metab).SNR(:);
        end
    end
    
    if strcmp(MRS_struct.p.Reference_compound,'H2O')
        out.water.FWHM = MRS_struct.out.(vox).water.FWHM(:);
        out.water.SNR  = MRS_struct.out.(vox).water.SNR(:);
    end
    out.NAA.FWHM   = MRS_struct.out.(vox).NAA.FWHM(:);
    out.NAA.SNR    = MRS_struct.out.(vox).NAA.SNR(:);
    
    if strcmpi(target{trg},'GABAGlx')
        if strcmp(MRS_struct.p.Reference_compound,'H2O')
            T = table(out.matlabVersion, out.GannetVersion, out.date, out.filename, round2(out.AvgDeltaF0), ...
                round2(out.GABA.ConcIU), round2(out.GABA.ConcCr), round2(out.GABA.FitError_W), round2(out.GABA.FitError_Cr), round2(out.GABA.FWHM), round2(out.GABA.SNR), ...
                round2(out.Glx.ConcIU), round2(out.Glx.ConcCr), round2(out.Glx.FitError_W), round2(out.Glx.FitError_Cr), round2(out.Glx.FWHM), round2(out.Glx.SNR), ...
                round2(out.water.FWHM), round2(out.water.SNR), round2(out.NAA.FWHM), round2(out.NAA.SNR), ...
                'VariableNames', {'MATLAB_version', 'Gannet_version', 'DateofAnalysis', 'Filename', 'Avg_Delta_F0', ...
                'GABA_ConcIU', 'GABA_ConcCr', 'GABA_FitError_W', 'GABA_FitError_Cr', 'GABA_FWHM','GABA_SNR', ...
                'Glx_ConcIU', 'Glx_ConcCr', 'Glx_FitError_W', 'Glx_FitError_Cr', 'Glx_FWHM','Glx_SNR', ...
                'H2O_FWHM', 'H2O_SNR', 'NAA_FWHM', 'NAA_SNR'});
        else
            T = table(out.matlabVersion, out.GannetVersion, out.date, out.filename, round2(out.AvgDeltaF0), ...
                round2(out.GABA.ConcCr), round2(out.GABA.FitError_Cr), round2(out.GABA.FWHM), round2(out.GABA.SNR), ...
                round2(out.Glx.ConcCr), round2(out.Glx.FitError_Cr), round2(out.Glx.FWHM), round2(out.Glx.SNR), ...
                round2(out.NAA.FWHM), round2(out.NAA.SNR), ...
                'VariableNames', {'MATLAB_version', 'Gannet_version', 'DateofAnalysis', 'Filename', 'Avg_Delta_F0', ...
                'GABA_ConcCr', 'GABA_FitError_Cr', 'GABA_FWHM','GABA_SNR', ...
                'Glx_ConcCr', 'Glx_FitError_Cr', 'Glx_FWHM','Glx_SNR', ...
                'NAA_FWHM', 'NAA_SNR'});
        end
    else
        if strcmp(MRS_struct.p.Reference_compound,'H2O')
            T = table(out.matlabVersion, out.GannetVersion, out.date, out.filename, round2(out.AvgDeltaF0), ...
                round2(out.(metab).ConcIU), round2(out.(metab).ConcCr), round2(out.(metab).FitError_W), round2(out.(metab).FitError_Cr), round2(out.(metab).FWHM), round2(out.(metab).SNR), ...
                round2(out.water.FWHM), round2(out.water.SNR), round2(out.NAA.FWHM), round2(out.NAA.SNR), ...
                'VariableNames', {'MATLAB_version', 'Gannet_version', 'DateofAnalysis', 'Filename', 'Avg_Delta_F0', ...
                [metab '_ConcIU'], [metab '_ConcCr'], [metab '_FitError_W'], [metab '_FitError_Cr'], [metab '_FWHM'], [metab '_SNR'], ...
                'H2O_FWHM', 'H2O_SNR', 'NAA_FWHM', 'NAA_SNR'});
        else
            T = table(out.matlabVersion, out.GannetVersion, out.date, out.filename, round2(out.AvgDeltaF0), ...
                round2(out.(metab).ConcCr), round2(out.(metab).FitError_Cr), round2(out.(metab).FWHM), round2(out.(metab).SNR), ...
                round2(out.NAA.FWHM), round2(out.NAA.SNR), ...
                'VariableNames', {'MATLAB_version', 'Gannet_version', 'DateofAnalysis', 'Filename', 'Avg_Delta_F0', ...
                [metab '_ConcCr'], [metab '_FitError_Cr'], [metab '_FWHM'], [metab '_SNR'], ...
                'NAA_FWHM', 'NAA_SNR'});
        end
    end
    
    
    % Break loop if function invoked in GannetFit
    if strcmpi(module,'fit')
        csv_name = ['MRS_struct_' target{trg} '_' vox '.csv'];
        writetable(T, csv_name);
        continue
    end
    
    %%% 2. Extract data from GannetSegment %%%
    
    if strcmp(MRS_struct.p.Reference_compound,'H2O')
        if strcmpi(target{trg},'GABAGlx')
            out.GABA.ConcIU_CSFcorr = MRS_struct.out.(vox).GABA.ConcIU_CSFcorr(:);
            out.Glx.ConcIU_CSFcorr = MRS_struct.out.(vox).Glx.ConcIU_CSFcorr(:);
        else
            out.(metab).ConcIU_CSFcorr = MRS_struct.out.(vox).(metab).ConcIU_CSFcorr(:);
        end
    end
    
    out.tissue.fGM  = MRS_struct.out.(vox).tissue.GMfra(:);
    out.tissue.fWM  = MRS_struct.out.(vox).tissue.WMfra(:);
    out.tissue.fCSF = MRS_struct.out.(vox).tissue.CSFfra(:);
    
    if strcmpi(target{trg},'GABAGlx')
        if strcmp(MRS_struct.p.Reference_compound,'H2O')
            U = table(round2(out.GABA.ConcIU_CSFcorr), round2(out.Glx.ConcIU_CSFcorr), ...
                round2(out.tissue.fGM), round2(out.tissue.fWM), round2(out.tissue.fCSF), ...
                'VariableNames', {'GABA_ConcIU_CSFcorr', 'Glx_ConcIU_CSFcorr', ...
                'fGM', 'fWM', 'fCSF'});
        else
            U = table(round2(out.tissue.fGM), round2(out.tissue.fWM), round2(out.tissue.fCSF), ...
                'VariableNames', {'fGM', 'fWM', 'fCSF'});
        end
    else
        if strcmp(MRS_struct.p.Reference_compound,'H2O')
            U = table(round2(out.(metab).ConcIU_CSFcorr), ...
                round2(out.tissue.fGM), round2(out.tissue.fWM), round2(out.tissue.fCSF), ...
                'VariableNames', {[metab '_ConcIU_CSFcorr'], ...
                'fGM', 'fWM', 'fCSF'});
        else
            U = table(round2(out.tissue.fGM), round2(out.tissue.fWM), round2(out.tissue.fCSF), ...
                'VariableNames', {'fGM', 'fWM', 'fCSF'});
        end
    end
    
    T = [T U]; %#ok<*AGROW>
    
    % Break loop if function invoked in GannetSegment
    if strcmpi(module,'segment')
        csv_name = ['MRS_struct_' target{trg} '_' vox '.csv'];
        writetable(T, csv_name);
        continue
    end
    
    %%% 3. Extract data from GannetQuantify %%%
    
    if strcmpi(target{trg},'GABAGlx')
        out.GABA.ConcIU_TissCorr              = MRS_struct.out.(vox).GABA.ConcIU_TissCorr(:);
        out.GABA.ConcIU_AlphaTissCorr         = MRS_struct.out.(vox).GABA.ConcIU_AlphaTissCorr(:);
        out.GABA.ConcIU_AlphaTissCorr_GrpNorm = MRS_struct.out.(vox).GABA.ConcIU_AlphaTissCorr_GrpNorm(:);
        out.Glx.ConcIU_TissCorr               = MRS_struct.out.(vox).Glx.ConcIU_TissCorr(:);
        out.Glx.ConcIU_AlphaTissCorr          = MRS_struct.out.(vox).Glx.ConcIU_AlphaTissCorr(:);
        out.Glx.ConcIU_AlphaTissCorr_GrpNorm  = MRS_struct.out.(vox).Glx.ConcIU_AlphaTissCorr_GrpNorm(:);
    else
        out.(metab).ConcIU_TissCorr              = MRS_struct.out.(vox).(metab).ConcIU_TissCorr(:);
        out.(metab).ConcIU_AlphaTissCorr         = MRS_struct.out.(vox).(metab).ConcIU_AlphaTissCorr(:);
        out.(metab).ConcIU_AlphaTissCorr_GrpNorm = MRS_struct.out.(vox).(metab).ConcIU_AlphaTissCorr_GrpNorm(:);
    end
    
    if strcmpi(target{trg},'GABAGlx')
        V = table(round2(out.GABA.ConcIU_TissCorr), round2(out.GABA.ConcIU_AlphaTissCorr), round2(out.GABA.ConcIU_AlphaTissCorr_GrpNorm), ...
            round2(out.Glx.ConcIU_TissCorr), round2(out.Glx.ConcIU_AlphaTissCorr), round2(out.Glx.ConcIU_AlphaTissCorr_GrpNorm), ...
            'VariableNames', {'GABA_ConcIU_TissCorr', 'GABA_ConcIU_AlphaTissCorr', 'GABA_ConcIU_AlphaTissCorr_GrpNorm', ...
            'Glx_ConcIU_TissCorr', 'Glx_ConcIU_AlphaTissCorr', 'Glx_ConcIU_AlphaTissCorr_GrpNorm'});
    else
        V = table(round2(out.(metab).ConcIU_TissCorr), round2(out.(metab).ConcIU_AlphaTissCorr), round2(out.(metab).ConcIU_AlphaTissCorr_GrpNorm), ...
            'VariableNames', {[metab '_ConcIU_TissCorr'], [metab '_ConcIU_AlphaTissCorr'], [metab '_ConcIU_AlphaTissCorr_GrpNorm']});
    end
    
    T = [T V];
    
    csv_name = ['MRS_struct_' target{trg} '_' vox '.csv'];
    writetable(T, csv_name);
    
end



