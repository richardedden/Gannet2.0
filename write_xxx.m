function new = write_xxx(string,old)
% Dynamic changing of regular expressions with regexprep.
% Overwrites strings within a string with the same number of X.
%
% Author: Dr. Georg Oeltzschner
%         Johns Hopkins University, 08/23/2016

% Strip unneeded " and white space.
old = regexprep(old, '\s*', ' ');
old = strrep(old,'"','');
old = strrep(old,' ','');
%
% Do actual overwriting.
new = strrep(string,old,repmat('X',[1 length(old)]));
end