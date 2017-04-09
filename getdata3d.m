function [ TrialbyVoxel ] = getdata3d(imagenames,varargin)
%[TrialByVoxel{mask}] = GETDATA3D(imagenames,maskinds)
%   getdata3d takes list of image names and mask indices and returns a
%   Trial x Voxel cell array
%   INPUT:
%   imagenames: can be a cell array of chars for multiple images
%       OR a single char for a single image
%   maskinds (optional) can be a cell array of inds for multiple
%       masks OR a numeric for a single mask

% check that input is correct
if isa(imagenames,'char')
    imgname = imagenames;
    clearvars imagenames
    imagenames{1}=imgname;
end

if ~isempty(varargin)
    
    if isa(varargin{1},'cell')
        maskinds = varargin{1};
    elseif isa(varargin{1},'numeric')
        maskinds{1} = varargin{1};
    end
else
    
    hdr = spm_vol(imagenames{1});
    maskinds{1} = (1:prod(hdr.dim))';
end

for mas = 1:length(maskinds)
    TrialbyVoxel{mas}=zeros(length(imagenames),length(maskinds{mas}));
end
for trial = 1:length(imagenames)
    hdr = spm_vol(imagenames{trial});
    img = spm_read_vols(hdr);
    for mas = 1:length(maskinds)
        TrialbyVoxel{mas}(trial,:) = img(maskinds{mas})';
%         S=whos('TrialbyVoxel');
%         if S.bytes/(10^9)>mem
%             error(['Call Exceeded Default Memory Threshold = ',num2str(mem),' gigs, If you have enough memory to proceed change memory by calling this function with ''memory''',',','numgigs where numgigs is a new memory threshold'])
%         end
    end
end

end

