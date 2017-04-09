function [ output_args ] = applywarp(warpprefix,moving,imagenameout,reference)
%APPLYWARP calls ANTS to warp your images
% applywarp(warpprefix,moving,imagenameout,reference)
%
% warpprefix - char with no image type extension
% moving - char with image type extension
% imagenameout - char with image type extension
% reference - char with image type extension
% 
% Don't forget to use absolute vs relative path when necessary and that
% ANTS doesn't like directories that use the ~ shortcut ex: ~/Documents/...

if ~isa(warpprefix,'char')
    error('All inputs must be chars');
end
if ~isa(moving,'char')
    error('All inputs must be chars');
end
if ~isa(imagenameout,'char')
    error('All inputs must be chars');
end
if ~isa(reference,'char')
    error('All inputs must be chars');
end

% if ~isempty(strfind(warpprefix,'.'))
%     error('No extensions in warp prefix')
% end
if ~ispc
callstr = '!/Users/colin/Applications/MyCode/antsbin/bin/WarpImageMultiTransform 3 ';
else
    callstr = '!WarpImageMultiTransform 3 ';
end
callstr=[callstr,moving,' '];
callstr=[callstr,imagenameout,' '];
callstr=[callstr,'-R ',reference,' --use-NN '];
warpimage = [warpprefix,'Warp.nii'];
warpaffine = [warpprefix,'Affine.txt'];
callstr=[callstr,warpimage,' '];
callstr=[callstr,warpaffine,' ']
eval(callstr)

end
