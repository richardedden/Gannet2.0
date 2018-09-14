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
                    if any(strcmpi(MRS_struct.p.target, {'GABA','Glx','GABAGlx'})) && strcmpi(MRS_struct.p.target2, 'GSH')
                        % 1=ExpA, 2=ExpB, 3=ExpC, 4=ExpD (MM: 171120)
                        MRS_struct.fids.ON_OFF  = repmat([1 1 0 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                        MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                    elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                        MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                    end
                    
                case 'Philips'
                    if any(strcmpi(MRS_struct.p.target, {'GABA','Glx','GABAGlx'})) && strcmpi(MRS_struct.p.target2, 'GSH')
                        % 1=?, 2=?, 3=?, 4=? (MM: 170703)
                        MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                        MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                    elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                        MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                    end
                    
                case 'Siemens_twix'
                    if any(strcmpi(MRS_struct.p.target, {'GABA','Glx','GABAGlx'})) && strcmpi(MRS_struct.p.target2, 'GSH')
                        % 1=?, 2=?, 3=?, 4=? (MM: 170703)
                        MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                        MRS_struct.fids.ON_OFF2 = repmat([1 1 0 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                    elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                        % This has not been tested with universal sequence -- 03142018 MGSaleh
                        % MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        % MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                    elseif strcmpi(MRS_struct.p.target, 'EtOH') && strcmpi(MRS_struct.p.target2, 'GSH')
                        MRS_struct.fids.ON_OFF  = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % EtOH
                        MRS_struct.fids.ON_OFF2 = repmat([0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                    end
                    
            end
            
        else % MEGA-PRESS
            
            if strcmpi(MRS_struct.p.vendor,'Philips') && strcmpi(MRS_struct.p.seqorig,'Philips')
                MRS_struct.fids.ON_OFF = [ones(1,size(MRS_struct.fids.data,2)/2) zeros(1,size(MRS_struct.fids.data,2)/2)];
            elseif strcmpi(MRS_struct.p.vendor,'Philips_data') % Hardcode for now. This repmat depends on the value entered for NSA
                if ceil(MRS_struct.p.LarmorFreq(MRS_struct.ii)) > 290
                    MRS_struct.fids.ON_OFF = repmat([ones(1,(MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows)) zeros(1,(MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows))],[1 size(MRS_struct.fids.data,2)/((MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows)*2)]); % GABA @ 7T % Changed by MGSaleh  -- 2017
                else
                    MRS_struct.fids.ON_OFF = repmat([1 1 0 0], [1 size(MRS_struct.fids.data,2)/4]);
                    %MRS_struct.fids.ON_OFF = repmat([1 0], [1 size(MRS_struct.fids.data,2)/2]); %This seems to work with the MR1 Philips_data
                end
            else
                MRS_struct.fids.ON_OFF = repmat([1 0], [1 size(MRS_struct.fids.data,2)/2]);
            end
            
        end
        
    case 'offfirst'
        
        if MRS_struct.p.HERMES % HERMES: GABAGlx or Lac and GSH (MM: 171120)
            
            switch MRS_struct.p.vendor
                
                case 'GE'
                    if any(strcmpi(MRS_struct.p.target, {'GABA','Glx','GABAGlx'})) && strcmpi(MRS_struct.p.target2, 'GSH')
                        % 1=?, 2=?, 3=?, 4=? (MM: 171120)
                        MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                        MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                    elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                        MRS_struct.fids.ON_OFF  = repmat([1 0 0 1], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                    end
                    
                case 'Philips'
                    if any(strcmpi(MRS_struct.p.target, {'GABA','Glx','GABAGlx'})) && strcmpi(MRS_struct.p.target2, 'GSH')
                        % 1=ExpC, 2=ExpB, 3=ExpA, 4=ExpD (MM: 170703)
                        MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                        MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                    elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                        MRS_struct.fids.ON_OFF  = repmat([1 0 0 1], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        MRS_struct.fids.ON_OFF2 = repmat([1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                    end
                    
                case 'Siemens_twix'
                    if any(strcmpi(MRS_struct.p.target, {'GABA','Glx','GABAGlx'})) && strcmpi(MRS_struct.p.target2, 'GSH')
                        % 1=?, 2=?, 3=?, 4=? (MM: 170703)
                        MRS_struct.fids.ON_OFF  = repmat([1 0 0 1], [1 size(MRS_struct.fids.data,2)/4]); % GABA
                        MRS_struct.fids.ON_OFF2 = repmat([0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                    elseif strcmpi(MRS_struct.p.target, 'GSH') && strcmpi(MRS_struct.p.target2, 'Lac')
                        % This has not been tested with universal sequence -- 03142018 MGSaleh
                        % MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                        % MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                    end
            end
            
        else % MEGA-PRESS
            
            if strcmpi(MRS_struct.p.vendor,'Philips') && strcmpi(MRS_struct.p.seqorig,'Philips')
                MRS_struct.fids.ON_OFF = [zeros(1,size(MRS_struct.fids.data,2)/2) ones(1,size(MRS_struct.fids.data,2)/2)];
            elseif strcmpi(MRS_struct.p.vendor,'Philips_data') % Hardcode for now. This repmat depends on the value entered for NSA
                if ceil(MRS_struct.p.LarmorFreq(MRS_struct.ii)) > 290
                    MRS_struct.fids.ON_OFF = repmat([zeros(1,(MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows)) ones(1,(MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows))],[1 size(MRS_struct.fids.data,2)/((MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows)*2)]); % GABA @ 7T % Changed by MGSaleh  -- 2017
                else
                    MRS_struct.fids.ON_OFF = repmat([0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
                    %                 MRS_struct.fids.ON_OFF = repmat([0 1], [1
                    %                 size(MRS_struct.fids.data,2)/2]); %This seems to work
                    %                 with the MR1 Philips_data
                end
            else
                MRS_struct.fids.ON_OFF = repmat([0 1], [1 size(MRS_struct.fids.data,2)/2]);
            end
            
        end
        
end

end

