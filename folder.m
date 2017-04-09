function [Out,Out2] = folder(wildcards,varargin)
% [Out] = folder(wildcards,filters)
%   Input a path/string of wildcards get folders
%   matching the wildcards 
% filters are function handles that work iteratively
% example
% [Out] = folder('*',@(S)isstrprop(S,'digit'),@any)
% returns any folder with a digit in it
% [Out] = folder('*',@(S)isstrprop(S,'digit'),@all)
% returns folders named with only digits
[wildpth] = fileparts(wildcards);
if ispc
    slash = '\';
else
    slash = '/';
end
x = dir(fullfile(wildcards));
if numel(x)==0
    Out={};
    return;
end
x = x([x(:).isdir]);
Out = {x(:).name}';
if ~isempty(varargin)
    [path,name,ext] = cellfun(@fileparts,Out,'UniformOutput',false);
    for i = 1:numel(varargin)
        name = cellfun(varargin{i},name,'UniformOutput',false);
    end
    Out = Out(cell2mat(name));
end
Out = cellfun(@(S)strcat(S,slash),Out,'uniformoutput',false);
Out = Out(~strcmp(Out(:), '../') & ~strcmp(Out(:), './'), :);
%Out2 = Out;
% if ispc
%     Out = cellfun(@(S)strcat(pwd,'\',S),Out,'uniformoutput',false);
% else
%     Out = cellfun(@(S)strcat(pwd,'/',S),Out,'uniformoutput',false);
% end
%Out = cellfun(@(S)fullfile(pwd,wildpth,S),Out,'uniformoutput',false);
end



