function params = getopt(optarglist, uplowcase, argin)

% GETOPT takes a list of all allowed parameters (plus optional values) as
% well as the function input parameters in argin. 
%
% IN:
%    optarglist .... Defaults for the optional arguments in a function
%                    given in a cell-array of the form
%                    {'param1', optional_value;
%                     'param2', []; % no ptional value}
%        
%    uplowcase ..... true: distinguish between upper/lower case letters
%                    false: do not distinguish between upper/lower case letters
%    argin ......... User provided input parameters. If accidentially one 
%                    parameter is specified more than once, the LAST 
%                    occurence will be taken as the value.
% OUT:
%    params ........ struct of parameters

% -------------------------------------------------------------------------
% project: uberall 
% -------------------------------------------------------------------------
% authors:
%    Matthias ROTH (MR), GIS, Uni Stuttgart 
%    Balaji Devaraju (BD), IfE, LUH, Hannover
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2018-11-27: MA, new parameter: 'display_warning'
%                    (more warnings on the screen for external users
%                     can be suppressed if your sure enough)
%    2018-07-25: MA, consider empty field for 'params.sub_wgs84'
%    2017-12-14: MA, checking 'params.sub_wgs84' for true/false
%    2017-06-06: BD, if case insensitive options then lower case optarglist
%    2016-02-25: BD, brush up help text
%    2015-02-17: MR, brush up help text
%    2014-10-22: MR, initial version
% -------------------------------------------------------------------------
% license:
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the  Free  Software  Foundation; either version 3 of the License, or
%    (at your option) any later version.
%  
%    This  program is distributed in the hope that it will be useful, but 
%    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
%    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
%    GNU General Public License for more details.
%  
%    You  should  have  received a copy of the GNU General Public License
%    along with Octave; see the file COPYING.  
%    If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

%% set predefined optional parameters
params = [];
for i = 1:size(optarglist, 1)
    if ischar(optarglist{i, 1})
        if uplowcase 
            param = optarglist{i, 1};
        else
            param = lower(optarglist{i, 1});
            optarglist{i, 1} = lower(optarglist{i, 1});
        end
        
        params = setfield(params, param, optarglist{i, 2});
    else
        error('Specify only strings as parameter names!');
    end
end
        
%% check for specified parameters    
if rem(length(argin) / 2, 1) > 0
    error('There is something wrong with your parameter list. Probably missing a parameter value.');
end
count = zeros(size(optarglist, 1), 1); % variable to check if parameters are specified only once

for i=1:2:length(argin)
    if uplowcase 
        param = argin{i};
    else
        param = lower(argin{i});
    end
    
    t = ismember(optarglist(:, 1), param);
    count = count + t;
    if sum(count > 1) > 0
        error(['Parameter "', param, '" specified more than once in the parameter list!']);
    end
    
    if sum(t) == 1
        params = setfield(params, param, argin{i + 1});
    elseif sum(t) > 1
        error(['Parameter "', param, '" specified more than once in the default parameter list!']);
    else
        error(['Unknown parameter: "', param, '"']);
    end
end


if isfield(params,'display_warning') == 0
    params.display_warning = true;
end
if islogical(params.display_warning)== false
    params.display_warning = true;
end



if isfield(params,'sub_wgs84') == 0
    params.sub_wgs84 = false
end
if islogical(params.sub_wgs84)== false
    switch params.sub_wgs84 
    case 1
        params.sub_wgs84 = true;
    case 0
         params.sub_wgs84 = false;
    otherwise
        error('please insert ''sub_wgs84'' as logical variable')
    end
end

