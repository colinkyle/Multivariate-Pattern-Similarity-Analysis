function [searchMaskInds,searchCenterInds,maskvol] = getsearchCenterInds(siz,maskinds,search_shape,expand)
%[searchvolinds,searchCenterInds,maskvol] = GETsearchCenterInds(siz,maskinds,search_shape,expand)
%
%   getsearchCenterInds computes the searchspace for a searchlight using ROI masks.
%   getsearchCenterInds will fill in any gaps around the masks to create one
%   solid volume.  Use expand to add padding arround that volume, I
%   recommend setting expand to be roughly the radius of your searchlight
%   (in voxels).
%   INPUTS:
%       siz - the size() of your functional images (hdr.dim if using spm_vol to load images)
%       maskinds - main output of getmasks() can be a cell array or plain array
%           of voxel indices corresponding to your search space
%       search_shape - a 3D array defining the shape of your searchlight.
%           Any non-zero values will be counted in the searchspace.
%           Start with something like this:
%           ex n=linspace(-1,1,5);
%              [x,y,z]=ndgrid(n,n,n);
%              search_sphere=(sqrt(x.^2 + y.^2 + z.^2)<=1);
%       expand - an integer used to pad the shape defined by maskinds
%   OUTPUTS:
%       searchMaskInds - searchmaskinds is a cell array of indices
%           corresponding to each searchlight location throught the
%           searchspace.  Formatting is similar to getmasks.
%       searchCenterInds - a vector with just the center ind of each search
%           volume.
%       maskvol - a volume with ones for each searchCenterInd

if ~isa(siz,'numeric')||numel(siz)~=3
    error('siz must be a vector of numel(3)')
end
if (~isa(search_shape,'numeric')&&~isa(search_shape,'logical'))||numel(size(search_shape))~=3||~all(mod(size(search_shape),2))
    error('search_shape must be a 3D matrix of ones and zeros with odd size for all dimensions see: help getsearchCenterInds')
end
if (rem(expand,1) ~= 0)||numel(expand)~=1
    error('expand must be an integer')
end

%% start function
maskvol = zeros(siz);
if isa(maskinds,'cell')
    if max(cellfun(@max,maskinds))>prod(siz)
        error('larger maskinds than in zeros(siz) please check that maskinds and siz are from same size images')
    end
    maskvol(cat(1,maskinds{:})) = 1;
elseif isa(maskinds,'numeric')
    if max(maskinds)>prod(siz)
        error('larger maskinds than in zeros(siz) please check that maskinds and siz are from same size images')
    end
    maskvol(maskinds(:)) = 1;
else
    error('maskinds must be a array or array of indices')
end

% fill the swiss cheese
%do x
for y = 1:siz(2)
    for z = 1:siz(3)
        x = maskvol(:,y,z);
        inds = find(x>0);
        if length(inds)>1
            inds = max((inds(1)-expand),1):min((inds(end)+expand),siz(1));
            x(inds) = 1;
            maskvol(:,y,z) = x;
            clearvars inds
        end
    end
end
clearvars x y z inds
%do y
for x = 1:siz(1)
    for z = 1:siz(3)
        y = maskvol(x,:,z);
        inds = find(y>0);
        if length(inds)>1
            inds = (inds(1)-expand):(inds(end)+expand);
            y(inds) = 1;
            maskvol(x,:,z) = y;
            clearvars inds
        end
    end
end
clearvars x y z inds
%do z
for x = 1:siz(1)
    for y = 1:siz(2)
        z = maskvol(x,y,:);
        inds = find(z>0);
        if length(inds)>1
            inds = (inds(1)-1):(inds(end)+1);
            inds(inds<1)=[];
            inds(inds>36)=[];
            z(inds) = 1;
            maskvol(x,y,:) = z;
            clearvars inds
        end
    end
end
clearvars x y z inds
% analyze search shape
search_shape(find(search_shape))=1;
searchcount = numel(find(search_shape));
% this convolution is how we tell where the searchlights can fit
if nargout<3
    maskvol = convn(maskvol,search_shape,'same');
    searchCenterInds = find(maskvol==searchcount);%all search centers
else
    maskvol1 = convn(maskvol,search_shape,'same');
    searchCenterInds = find(maskvol1==searchcount);%all search centers
end
%for each searchcenter we want to return the entire searchlight volume
[I,J,K] = ind2sub(size(search_shape),find(search_shape));
xyz = [I,J,K];
xyz = bsxfun(@minus,xyz,xyz(round(numel(I)/2),:));%strange things might happen if the center of your search_shape is not a one
[i,j,k]=ind2sub(siz,searchCenterInds);
%searchmaskinds = zeros(searchcount,numel(searchCenterInds));
for s = 1:numel(searchCenterInds)
    XYZ = bsxfun(@plus,xyz,[i(s),j(s),k(s)]);%finds coordinates relative to each search center
    searchMaskInds{s} = XYZ(:,1)+siz(1)*(XYZ(:,2)-1)+siz(1)*siz(2)*(XYZ(:,3)-1);%convert coordinates to indices
end

end

