function [] = callspmsegmentation(anatimage)




%get file handling for spm sorted out
[spmhome spmdir spmect] = fileparts(which('spm'));

matlabbatch{1}.spm.spatial.preproc.data = {[anatimage ',1']};
matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{1}.spm.spatial.preproc.output.cleanup = 0;
matlabbatch{1}.spm.spatial.preproc.opts.tpm = {
                                               [spmhome '/tpm/grey.nii']
                                               [spmhome '/tpm/white.nii']
                                               [spmhome '/tpm/csf.nii']
                                               };
matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2
                                                 2
                                                 2
                                                 4];
matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};
matlabbatch{2}.spm.tools.preproc8.channel.vols = {[anatimage ',1']};
matlabbatch{2}.spm.tools.preproc8.channel.biasreg = 0.0001;
matlabbatch{2}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{2}.spm.tools.preproc8.channel.write = [0 0];
matlabbatch{2}.spm.tools.preproc8.tissue(1).tpm = {[spmhome '/toolbox/Seg/TPM.nii,1']};
matlabbatch{2}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{2}.spm.tools.preproc8.tissue(1).native = [1 0];
matlabbatch{2}.spm.tools.preproc8.tissue(1).warped = [0 0];
matlabbatch{2}.spm.tools.preproc8.tissue(2).tpm = {[spmhome '/toolbox/Seg/TPM.nii,2']};
matlabbatch{2}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{2}.spm.tools.preproc8.tissue(2).native = [1 0];
matlabbatch{2}.spm.tools.preproc8.tissue(2).warped = [0 0];
matlabbatch{2}.spm.tools.preproc8.tissue(3).tpm = {[spmhome '/toolbox/Seg/TPM.nii,3']};
matlabbatch{2}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{2}.spm.tools.preproc8.tissue(3).native = [1 0];
matlabbatch{2}.spm.tools.preproc8.tissue(3).warped = [0 0];
matlabbatch{2}.spm.tools.preproc8.tissue(4).tpm = {[spmhome '/toolbox/Seg/TPM.nii,4']};
matlabbatch{2}.spm.tools.preproc8.tissue(4).ngaus = 3;
matlabbatch{2}.spm.tools.preproc8.tissue(4).native = [1 0];
matlabbatch{2}.spm.tools.preproc8.tissue(4).warped = [0 0];
matlabbatch{2}.spm.tools.preproc8.tissue(5).tpm = {[spmhome '/toolbox/Seg/TPM.nii,5']};
matlabbatch{2}.spm.tools.preproc8.tissue(5).ngaus = 4;
matlabbatch{2}.spm.tools.preproc8.tissue(5).native = [1 0];
matlabbatch{2}.spm.tools.preproc8.tissue(5).warped = [0 0];
matlabbatch{2}.spm.tools.preproc8.tissue(6).tpm = {[spmhome '/toolbox/Seg/TPM.nii,6']};
matlabbatch{2}.spm.tools.preproc8.tissue(6).ngaus = 2;
matlabbatch{2}.spm.tools.preproc8.tissue(6).native = [0 0];
matlabbatch{2}.spm.tools.preproc8.tissue(6).warped = [0 0];
matlabbatch{2}.spm.tools.preproc8.warp.mrf = 0;
matlabbatch{2}.spm.tools.preproc8.warp.reg = 4;
matlabbatch{2}.spm.tools.preproc8.warp.affreg = 'mni';
matlabbatch{2}.spm.tools.preproc8.warp.samp = 3;
matlabbatch{2}.spm.tools.preproc8.warp.write = [0 0];

spm_jobman('initcfg')
spm_jobman('run',matlabbatch)

end

