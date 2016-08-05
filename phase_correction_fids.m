function out = phase_correction_fids(data,water_data)

% Written and updated by MGSaleh, University of Cape Town, 2015.
% The function subtracts the phase of the unsuppressed (water_data) data
% from the phase of the suppressed (data) data.

% in=data;
K=abs(water_data(1,:));  % Only takes the first coloumn of data if water data are not averaged  -- MGSaleh
Kphase=phase(water_data(1,:));

sz_mat=size(data);
fids=zeros(sz_mat);
% specs=zeros(sz_mat);



for kk=1:sz_mat(1)
%     for ll=1:sz_mat(2) %(sz_mat(2)*sz_mat(3))
% figure,plot(real(fftshift(fft(data(kk,:)))))
% xlim([900 1500])
% hold on


    K_data=abs(data(kk,:));
    Kphase_data=phase(data(kk,:));
    
    Kphase_corr=(Kphase_data)-Kphase;
    
    fids(kk,:)=K_data.* exp(1i*Kphase_corr);
    %         specs(:,ll)=fftshift(ifft(fids(:,ll,kk),[],data.dims.t),data.dims.t);
    
%     end
% figure,plot(real(fftshift(fft(fids(kk,:)))))
% xlim([900 1500])


end



out = fids;


% end


% close all
%
% figure,plot(in.ppm,(real(((fftshift(ifft((in.fids(:,:,1)),[],in.dims.t),in.dims.t))))))
% xlim([1.5 4])
% figure,plot(in.ppm,(real(((fftshift(ifft((fids(:,:,1)),[],in.dims.t),in.dims.t))))))
% xlim([1.5 4])
% figure,plot(in.ppm,(real(((fftshift(ifft((in.fids(:,:,2)),[],in.dims.t),in.dims.t))))))
% xlim([1.5 4])
% figure,plot(in.ppm,(real(((fftshift(ifft((fids(:,:,2)),[],in.dims.t),in.dims.t))))))
% xlim([1.5 4])
%
% pause
% close all;
%
% out.sz=in.sz;
% out.ppm=in.ppm;
% out.t=in.t;
% out.spectralwidth=in.spectralwidth; %Changed by MSaleh
% out.dwelltime=in.dwelltime; %Defined differently by MSaleh
% out.txfrq=in.txfrq;  %txfrq; %Changed by MSaleh
% % out.date=date; %Changed by MSaleh
% out.dims=in.dims;
% out.Bo=in.Bo;
% out.averages=in.averages;
% out.rawAverages=in.rawAverages;
% out.subspecs=in.subspecs;
% out.rawSubspecs=in.rawSubspecs;
% clear out.fids out.specs;
% out.fids=fids;
% out.specs=specs;
%
% %FILLING IN THE FLAGS
% out.flags=in.flags;
% out.flags.eddy_corrected=1;
% end
%
% % if in_w.eddy_correction
% %             disp('Perfoming eddy correction')
% %         size(out1_ls.fids)
% %         size(in_w.fids)
% %                 out1_ls.fids(:,1)=klose_eddy_correction(out1_ls.fids(:,1),in_w.fids(:,1,2));
% %                 out1_ls.fids(:,2)=klose_eddy_correction(out1_ls.fids(:,2),in_w.fids(:,1,2));
% %          end