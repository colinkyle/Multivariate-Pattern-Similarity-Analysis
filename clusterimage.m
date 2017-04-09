function [ClusterInds, MaxSize, NumSigVox ] = clusterimage(img,thresh)
% [ClusterInds, MaxSize, MaxCount ] = CLUSTERIMAGE(IMG,THRESH)
%    clusters your data above a threshold and returns the
%    resulting indices(sorted by size), the max cluster size size, and the 
%    number of significant voxels
%    IMG - can either be a char (nifti file name), or a 3d image array

if isa(img,'char')
    hdr = spm_vol(img);
    clearvars img
    img = spm_read_vols(hdr);
end

img = img>=thresh;
CC = bwconncomp(img);
AllSize = cellfun(@length,CC.PixelIdxList);
[~,I] = sort(AllSize,'descend');
ClusterInds = CC.PixelIdxList(I);
MaxSize = max(cellfun(@length,CC.PixelIdxList));
NumSigVox = sum(cellfun(@length,CC.PixelIdxList));
if isempty(MaxSize)
    MaxSize = 0;
end

end

