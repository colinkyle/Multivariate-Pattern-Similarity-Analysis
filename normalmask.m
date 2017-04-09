function [success] = normalmask(filenameout,maskinds,basevals,hdr)
%NORMALMASK takes an output file name, maskinds output from getmasks, and a vector of base
%values the same length as maskinds and an image header from spm_vol
% It writes an image named filenameout containing normal distributions for
% each mask and a larger distribution from the center of all masks through
% the rest of the image.  basevals should all be greater than 100 and the
% same for all subjects.
%
%ex normalmask(filenameout,maskinds,basevals,hdr)
%   filenameout: 'string.nii'
%   maskinds: {1xN_masks cell}
%   basevals: [1xN_masks double]
%   hdr: image header file (from spm_vol)

%   Detailed explanation goes here
%%
if ~isa(filenameout,'char')
    error('output file name must be char');
end
if ~isa(maskinds,'cell')&&~isa(maskinds,'numeric')
    error('maskinds must be numeric or cell of numerics')
end
if ~isa(basevals,'numeric') || (numel(basevals) ~= numel(maskinds))
    error('basevals must be a numeric with length=length(maskinds)')
end

img = zeros(hdr.dim);
for mas = 1:length(maskinds)
    tempimg = zeros(hdr.dim);
    tempimg(maskinds{mas})=1;
    [I1,J1,K1] = ind2sub(hdr.dim,maskinds{mas});
    mins= min([I1,J1,K1]);
    maxs = max([I1,J1,K1]);
    %largest spread of X,Y,Z becomes basis for distrubution
    maxdim = max(maxs-mins);
    %Get Centroid
    stats = regionprops(tempimg,'centroid');
    %Create and normalize distribution (distribution span set to 5)
    x = [0:max(hdr.dim)];
    normdist = normpdf(x,0,round(maxdim/7));
    normdist = 5*normdist/max(normdist);
    %Go through inds set value of mask to
    %normdist(radius_from_centroid)+(constant to distinguish mask set)
    for i = 1:length(maskinds{mas})
        [I,J,K] = ind2sub(hdr.dim,maskinds{mas}(i));
        dist = round(pdist([[stats.Centroid([2 1 3])];[I,J,K]],'euclidean'));
        img(maskinds{mas}(i))=normdist(dist+1)+basevals(mas);
    end
end
inds = find(img==0);
[I1,J1,K1] = ind2sub(hdr.dim,find(img==0));
mins= min([I1,J1,K1]);
maxs = max([I1,J1,K1]);
%largest spread of X,Y,Z becomes basis for distrubution
maxdim = max(maxs-mins);
%Get Centroid
dummyimg =img;
dummyimg(img~=0)=1;
stats = regionprops(dummyimg,'centroid');
%Create and normalize distribution (distribution span set to 5)
x = [0:max(hdr.dim)];
normdist = normpdf(x,0,round(maxdim/7));
normdist = 5*normdist/max(normdist);
%Go through inds set value of mask to
%normdist(radius_from_centroid)+(constant to distinguish mask set)
for i = 1:length(inds)
    [I,J,K] = ind2sub(hdr.dim,inds(i));
    dist = round(pdist([[stats.Centroid([2 1 3])];[I,J,K]],'euclidean'));
    img(inds(i))=normdist(dist+1)+5;
end
hdr.fname = filenameout;
spm_write_vol(hdr,img);
end

