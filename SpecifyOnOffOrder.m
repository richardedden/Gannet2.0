function MRS_struct = SpecifyOnOffOrder(MRS_struct)
%   Determines the order of OFF and ON for HERMES and MEGA-PRESS.
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-02-24)
%       goeltzs1@jhmi.edu
%
%   History:
%       2018-02-24: First version.

% [1 = GABA_ON, 0 = GABA_OFF]
switch MRS_struct.p.ONOFForder
    
    case 'onfirst'
        
        if MRS_struct.p.HERMES % HERMES: GABAGlx or Lac and GSH (MM: 171120)
            
            switch MRS_struct.p.vendor
                
                case 'GE'
                    if strcmpi(MRS_struct.p.target, 'GABAGlx') && strcmpi(MRS_struct.p.target2, 'GSH')
                        % 1=ExpA, 2=ExpB, 3=ExpC, 4=ExpD (MM: 171120)
                        MRS_struct.fids.ON_OFF  = repmat([1 1 0 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                        MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                    elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                        MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                    end
                    
                case 'Philips'
                    if strcmpi(MRS_struct.p.target, 'GABAGlx') && strcmpi(MRS_struct.p.target2, 'GSH')
                        % 1=?, 2=?, 3=?, 4=? (MM: 170703)
                        MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                        MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                    elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                        MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                    end
                
                case 'Siemens_twix'
                    if strcmpi(MRS_struct.p.target, 'GABAGlx') && strcmpi(MRS_struct.p.target2, 'GSH')
                        % 1=?, 2=?, 3=?, 4=? (MM: 170703)
                        MRS_struct.fids.ON_OFF  = repmat([1 1 0 0 0 0 1 1], [1 size(MRS_struct.fids.data,2)/8]); % GABA
                        MRS_struct.fids.ON_OFF2 = repmat([0 0 0 0 1 1 1 1], [1 size(MRS_struct.fids.data,2)/8]); % GSH
                    elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac') %This has not been tested with universal sequence -- 03142018 MGSaleh
                        % MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        % MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                    end
                    
            end
            
        else % MEGAPRESS
            
            if strcmpi(MRS_struct.p.vendor,'Philips') && strcmpi(MRS_struct.p.seqorig,'Philips')
                MRS_struct.fids.ON_OFF = [ones(1,size(MRS_struct.fids.data,2)/2) zeros(1,size(MRS_struct.fids.data,2)/2)];
            elseif strcmpi(MRS_struct.p.vendor,'Philips_data') % Hardcode for now. This repmat depends on the value entered for NSA
                MRS_struct.fids.ON_OFF = repmat([1 1 0 0], [1 size(MRS_struct.fids.data,2)/4]);
            elseif strcmpi(MRS_struct.p.vendor,'Siemens_twix') && strcmp(MRS_struct.p.seq_string,'mgs_svs_ed') %This is a condition to check whether this is a universal sequence -- 03162018  MGSaleh
                MRS_struct.fids.ON_OFF = repmat([1 1 0 0],[1 size(MRS_struct.fids.data,2)/4]);
            else
                MRS_struct.fids.ON_OFF = repmat([1 0], [1 size(MRS_struct.fids.data,2)/2]);
            end

        end
        
        if strcmp(MRS_struct.p.seq_string,'mgs_svs_ed') %This is a condition to check whether this is a universal sequence -- 03162018  MGSaleh
            %Zero-order correction based on Cr -- 03162018  MGSaleh
            MRS_struct = cr_phase_corr_univ_seq(MRS_struct);
        end

    case 'offfirst'
        
        if MRS_struct.p.HERMES % HERMES: GABAGlx or Lac and GSH (MM: 171120)
            
            switch MRS_struct.p.vendor
                
                case 'GE'
                    if strcmpi(MRS_struct.p.target, 'GABAGlx') && strcmpi(MRS_struct.p.target2, 'GSH')
                        % 1=?, 2=?, 3=?, 4=? (MM: 171120)
                        MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                        MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                    elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                        MRS_struct.fids.ON_OFF  = repmat([1 0 0 1], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                    end
                    
                case 'Philips'
                    if strcmpi(MRS_struct.p.target, 'GABAGlx') && strcmpi(MRS_struct.p.target2, 'GSH')
                        % 1=ExpC, 2=ExpB, 3=ExpA, 4=ExpD (MM: 170703)
                        MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                        MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                    elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                        MRS_struct.fids.ON_OFF  = repmat([1 0 0 1], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                    end
            end
            
        else % MEGAPRESS
            if strcmpi(MRS_struct.p.vendor,'Philips') && strcmpi(MRS_struct.p.seqorig,'Philips')
                MRS_struct.fids.ON_OFF = [zeros(1,size(MRS_struct.fids.data,2)/2) ones(1,size(MRS_struct.fids.data,2)/2)];
            elseif strcmpi(MRS_struct.p.vendor,'Philips_data') % Hardcode for now. This repmat depends on the value entered for NSA
                MRS_struct.fids.ON_OFF = repmat([0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
            else
                MRS_struct.fids.ON_OFF = repmat([0 1], [1 size(MRS_struct.fids.data,2)/2]);
            end
        end
        
end

end

