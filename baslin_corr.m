close all

z=abs(MRS_HERMES_GABA_GSH.spec.freq-4);
lowerbound=find(min(z)==z);
z=abs(MRS_HERMES_GABA_GSH.spec.freq-2.8);        %2.75
upperbound=find(min(z)==z);
freqbounds=lowerbound:upperbound;


ii=1;

yy=size(MRS_HERMES_GABA_GSH.spec.GSH.diff(ii,:));

low_val=min(real(MRS_HERMES_GABA_GSH.spec.GSH.diff(ii,[freqbounds])));

% MRS_HERMES_GABA_GSH.spec.GSH.diff_new(ii,:)= MRS_HERMES_GABA_GSH.spec.GSH.diff(ii,:); 
MRS_HERMES_GABA_GSH.spec.GSH.diff_nobas_corr(ii,:)=MRS_HERMES_GABA_GSH.spec.GSH.diff(ii,:);

MRS_HERMES_GABA_GSH.spec.GSH.diff(ii,:) = complex(real(MRS_HERMES_GABA_GSH.spec.GSH.diff(ii,:)) - 1*low_val*(ones(yy)), imag(MRS_HERMES_GABA_GSH.spec.GSH.diff(ii,:)));

figure,plot(MRS_HERMES_GABA_GSH.spec.freq,real(MRS_HERMES_GABA_GSH.spec.GSH.diff_nobas_corr(ii,:)))
% xlim([freqbounds])

figure,plot(MRS_HERMES_GABA_GSH.spec.freq,real(MRS_HERMES_GABA_GSH.spec.GSH.diff(ii,:)))
% xlim([freqbounds])