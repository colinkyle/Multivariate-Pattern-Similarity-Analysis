function [maskinds,dimensions,hdr] = getmasks(filename,varargin)
% [INDS,DIMENSION,hdr] = GETMASK(FILENAME,[MASKVAL])
%    returns indices of mask volume, dimensions of image, and the mask's
%    header (from spm_vol)
%    FILENAME - can be char or cell of chars
%       * Filename can also be a volume
%    MASKVAL - optional, default is find(mask_image>0) or
%      find(mask_image==MASKVAL)
%    *** This function depends on spm_vol and spm_read_vols (you need SPM!)
if isempty(which('spm_vol'))
    error('You need to install SPM or add it to your path!')
end
if ~isnumeric(filename)
if isa(filename,'cell')
    for ii = 1:length(filename)
        hdr = spm_vol(filename{ii});
        img = spm_read_vols(hdr);
        if length(varargin)
            maskinds{ii} = find(img==varargin{1});
        else
            maskinds{ii} = find(img>0);
        end
        dimensions{ii} = size(img);
    end
else
    hdr = spm_vol(filename);
    img = spm_read_vols(hdr);
    if length(varargin)
        maskinds = find(img==varargin{1});
    else
        maskinds = find(img>0);
    end
    dimensions = size(img);
end
else
    if length(varargin)
        for ii = 1:numel(varargin{1})
            maskinds{ii} = find(filename==varargin{1}(ii));
        end
    else
        maskinds = find(filename>0);
    end
end
end

