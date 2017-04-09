function [Out,Out2] = file(wildcards,varargin)
% [Out] = file(wildcards,filters)
%   Input a path/string of wildcards get filenames
%   matching the wildcards 
% filters are function handles that work iteratively
% example
% [Out] = file('*',@(S)isstrprop(S,'digit'),@any)
% returns any file with a digit in it
% [Out] = file('*',@(S)isstrprop(S,'digit'),@all)
% returns files named with only digits
[wildpth] = fileparts(wildcards);
x = dir(fullfile(wildcards));

x = x(~[x(:).isdir]);
Out = {x(:).name};
Out = Out';
Out = sort_nat(Out);
if ~isempty(varargin)
    [path,name,ext] = cellfun(@fileparts,Out,'UniformOutput',false);
    for i = 1:numel(varargin)
        name = cellfun(varargin{i},name,'UniformOutput',false);
    end
    Out = Out(cell2vec(name));
end
% Out2 = Out;
% 
% Out = cellfun(@(S)fullfile(pwd,wildpth,S),Out,'uniformoutput',false);

end


