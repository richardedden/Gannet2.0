function MRS_struct = cr_phase_corr_univ_seq(MRS_struct)

% Performing correction on the HERMES data using Cr phase of the 1st average of every experiment  --  03152018 MGSaleh
close all

% Create some new variables and use them to perform the intended
% correction. NOTE: The MRS_struct.fids.data gets overwritten at the end of
% this code -- MGSaleh
FullData = MRS_struct.fids.data;
% figure(1),plot(real(FullData))
ZeroFillTo(MRS_struct.ii) = round(32768/2000*MRS_struct.p.sw(MRS_struct.ii));
zf = ZeroFillTo(MRS_struct.ii)/MRS_struct.p.npoints(MRS_struct.ii);
time = (1:1:size(FullData,1))/MRS_struct.p.sw(MRS_struct.ii);
% Line-broadening, zero-filling and FFT
FullData2 = FullData .* repmat((exp(-time'*MRS_struct.p.LB*pi)), [1 size(FullData,2)]);
AllFramesFT_old = fftshift(fft(FullData2,ZeroFillTo(MRS_struct.ii),1),1);
% Work out frequency scale
freqrange = MRS_struct.p.sw(MRS_struct.ii)/MRS_struct.p.LarmorFreq(MRS_struct.ii);
freq_data = (ZeroFillTo(MRS_struct.ii)+1-(1:1:ZeroFillTo(MRS_struct.ii)))/ZeroFillTo(MRS_struct.ii)*freqrange+4.68-freqrange/2.0;

% Save the first average of every experiment. The first four variables should work for MEGA-PRESS as well -- MGSaleh
ONGABAOFFGSH1=abs(AllFramesFT_old(:,1));  % This is equivalent of ONGABA1 in MEGA-PRESS -- MGSaleh
ONGABAOFFGSH2=abs(AllFramesFT_old(:,2));  % This is equivalent of ONGABA2 in MEGA-PRESS -- MGSaleh
OFFGABAOFFGSH1=abs(AllFramesFT_old(:,3)); % This is equivalent of OFFGABA1 in MEGA-PRESS -- MGSaleh
OFFGABAOFFGSH2=abs(AllFramesFT_old(:,4)); % This is equivalent of OFFGABA2 in MEGA-PRESS -- MGSaleh
if MRS_struct.p.HERMES % We run through this code if HERMES -- MGSaleh
    OFFGABAONGSH1=abs(AllFramesFT_old(:,5));
    OFFGABAONGSH2=abs(AllFramesFT_old(:,6));
    ONGABAONGSH1=abs(AllFramesFT_old(:,7));
    ONGABAONGSH2=abs(AllFramesFT_old(:,8));
end

% Define a range for the search of peak position of Cr --  MGSaleh
z=abs(freq_data-3.1);
lowerbound=find(min(z)==z);
z=abs(freq_data-2.95);%2.75
upperbound=find(min(z)==z);
freqbounds=lowerbound:upperbound;

% Determine the positinon of the peak position of Cr --  MGSaleh
ONGABAOFFGSH1max  = max(ONGABAOFFGSH1(freqbounds,1));
ONGABAOFFGSH2max  = max(ONGABAOFFGSH2(freqbounds,1));
OFFGABAOFFGSH1max = max(OFFGABAOFFGSH1(freqbounds,1));
OFFGABAOFFGSH2max = max(OFFGABAOFFGSH2(freqbounds,1));
if MRS_struct.p.HERMES % We run through this code if HERMES -- MGSaleh
    OFFGABAONGSH1max  = max(OFFGABAONGSH1(freqbounds,1));
    OFFGABAONGSH2max  = max(OFFGABAONGSH2(freqbounds,1));
    ONGABAONGSH1max   = max(ONGABAONGSH1(freqbounds,1));
    ONGABAONGSH2max   = max(ONGABAONGSH2(freqbounds,1));
end

% Determine the complex value at the max position -- MGSaleh
ph1=find(ONGABAOFFGSH1==ONGABAOFFGSH1max);
ph2=find(ONGABAOFFGSH2==ONGABAOFFGSH2max);
ph3=find(OFFGABAOFFGSH1==OFFGABAOFFGSH1max);
ph4=find(OFFGABAOFFGSH2==OFFGABAOFFGSH2max);
if MRS_struct.p.HERMES % We run through this code if HERMES -- MGSaleh
    ph5=find(OFFGABAONGSH1==OFFGABAONGSH1max);
    ph6=find(OFFGABAONGSH2==OFFGABAONGSH2max);
    ph7=find(ONGABAONGSH1==ONGABAONGSH1max);
    ph8=find(ONGABAONGSH2==ONGABAONGSH2max);
end

% Determine the phase at the same position -- MGSaleh
phas_cr1 = phase(AllFramesFT_old(ph1,1));
phas_cr2 = phase(AllFramesFT_old(ph2,2));
phas_cr3 = phase(AllFramesFT_old(ph3,3));
phas_cr4 = phase(AllFramesFT_old(ph4,4));
if MRS_struct.p.HERMES % We run through this code if HERMES -- MGSaleh
    phas_cr5 = phase(AllFramesFT_old(ph5,5));
    phas_cr6 = phase(AllFramesFT_old(ph6,6));
    phas_cr7 = phase(AllFramesFT_old(ph7,7));
    phas_cr8 = phase(AllFramesFT_old(ph8,8));
end

% Visually determine if the corrections works -- MGSaleh
figure(2)
if MRS_struct.p.HERMES % We run through this code if HERMES -- MGSaleh
    n_plots = 8;
else
    n_plots = 4;
end

subplot(2,n_plots,1),plot(freq_data,real(fftshift(fft(FullData(:,1:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*0))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAON and GSHOFF before correction'),xlabel('ppm'),ylabel('a.u')
subplot(2,n_plots,1+n_plots),plot(freq_data,real(fftshift(fft(FullData(:,1:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*(phas_cr1)))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAON and GSHOFF after correction'),xlabel('ppm'),ylabel('a.u')
subplot(2,n_plots,2),plot(freq_data,real(fftshift(fft(FullData(:,2:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*0))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAON and GSHOFF before correction'),xlabel('ppm'),ylabel('a.u')
subplot(2,n_plots,2+n_plots),plot(freq_data,real(fftshift(fft(FullData(:,2:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*(phas_cr2)))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAON and GSHOFF after correction'),xlabel('ppm'),ylabel('a.u')
%
subplot(2,n_plots,3),plot(freq_data,real(fftshift(fft(FullData(:,3:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*0))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAOFF and GSHOFF before correction'),xlabel('ppm'),ylabel('a.u')
subplot(2,n_plots,3+n_plots),plot(freq_data,real(fftshift(fft(FullData(:,3:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*(phas_cr3)))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAOFF and GSHOFF after correction'),xlabel('ppm'),ylabel('a.u')
subplot(2,n_plots,4),plot(freq_data,real(fftshift(fft(FullData(:,4:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*0))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAOFF and GSHOFF before correction'),xlabel('ppm'),ylabel('a.u')
subplot(2,n_plots,4+n_plots),plot(freq_data,real(fftshift(fft(FullData(:,4:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*(phas_cr4)))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAOFF and GSHOFF after correction'),xlabel('ppm'),ylabel('a.u')

if MRS_struct.p.HERMES % We run through this code if HERMES -- MGSaleh
    subplot(2,n_plots,5),plot(freq_data,real(fftshift(fft(FullData(:,5:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*0))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAOFF and GSHON before correction'),xlabel('ppm'),ylabel('a.u')
    subplot(2,n_plots,5+n_plots),plot(freq_data,real(fftshift(fft(FullData(:,5:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*(phas_cr5)))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAOFF and GSHON after correction'),xlabel('ppm'),ylabel('a.u')
    subplot(2,n_plots,6),plot(freq_data,real(fftshift(fft(FullData(:,6:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*0))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAOFF and GSHON before correction'),xlabel('ppm'),ylabel('a.u')
    subplot(2,n_plots,6+n_plots),plot(freq_data,real(fftshift(fft(FullData(:,6:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*(phas_cr6)))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAOFF and GSHON after correction'),xlabel('ppm'),ylabel('a.u')
    
    subplot(2,n_plots,7),plot(freq_data,real(fftshift(fft(FullData(:,7:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*0))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAON and GSHON before correction'),xlabel('ppm'),ylabel('a.u')
    subplot(2,n_plots,7+n_plots),plot(freq_data,real(fftshift(fft(FullData(:,7:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*(phas_cr7)))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAON and GSHON after correction'),xlabel('ppm'),ylabel('a.u')
    subplot(2,n_plots,8),plot(freq_data,real(fftshift(fft(FullData(:,8:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*0))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAON and GSHON before correction'),xlabel('ppm'),ylabel('a.u')
    subplot(2,n_plots,8+n_plots),plot(freq_data,real(fftshift(fft(FullData(:,8:8:end),ZeroFillTo(MRS_struct.ii),1),1).*exp(-1i*(phas_cr8)))),set(gca,'Xdir','reverse'),xlim([1 4]),title('Edit GABAON and GSHON after correction'),xlabel('ppm'),ylabel('a.u')
end

% Apply the correction --  MGSaleh
FullData(:,1:n_plots:end)=FullData(:,1:n_plots:end).*exp(-1i*(phas_cr1));
FullData(:,2:n_plots:end)=FullData(:,2:n_plots:end).*exp(-1i*(phas_cr2));
FullData(:,3:n_plots:end)=FullData(:,3:n_plots:end).*exp(-1i*(phas_cr3));
FullData(:,4:n_plots:end)=FullData(:,4:n_plots:end).*exp(-1i*(phas_cr4));
if MRS_struct.p.HERMES % We run through this code if HERMES -- MGSaleh
    FullData(:,5:n_plots:end)=FullData(:,5:n_plots:end).*exp(-1i*(phas_cr5));
    FullData(:,6:n_plots:end)=FullData(:,6:n_plots:end).*exp(-1i*(phas_cr6));
    FullData(:,7:n_plots:end)=FullData(:,7:n_plots:end).*exp(-1i*(phas_cr7));
    FullData(:,8:n_plots:end)=FullData(:,8:n_plots:end).*exp(-1i*(phas_cr8));
end

% Save the data back to the original variable --  MGSaleh
MRS_struct.fids.data = FullData;

% figure(3),plot(freq_data,real(fftshift(fft(out.fids.data,ZeroFillTo(MRS_struct.ii),1),1))),set(gca,'Xdir','reverse'),xlim([1 4])
% figure(3),plot(real(FullData))

% disp('Paused. Press enter to continue')
% 
% pause

close all

end