function [SimilarityVec] = mpsvec(TrialByVoxel1,TrialByVoxel2,varargin)
%[SimilarityVec] = MPSVEC(TrialByVoxel1,TrialByVoxel2)
%MPSVEC Returns the r values defined by comparisons between matching
%   trials between two blocks of data. (i.e. Trial1a vs. Trial2a, Trial1b
%   vs. Trial2b, ...)  Therefore, TrialByVoxel1 and TrialByVoxel2 must have
%   the same number of trials.
%   
%   Optional input argument name-value pairs:
%      'outlier',num - removes outliers [num] standard deviations from voxel
%           mean.  Default behavior num=Inf.  Voxels that deviate by more than 
%           num standard deviations will be treated as NaNs and thus removed
%           from the analysis completely unless 'replaceNaN' is true.
%
%      'replaceNaN',bool - this sets the behavior for dealing with NaNs/outliers.
%           By default (false), this function removes voxels containing NaNs/outliers 
%           from both blocks of data.  This behavior is not always desireable. 
%           For instance, it can cause errors when sending a searchlight 
%           through noisy data.  Setting replaceNaN to true will cause NaNs/outliers 
%           to be replaced by the mean voxel value for a given trial.  This 
%           prohibits a voxel from contributing positively or negatively to 
%           a correlation. Thus, setting 'replaceNaN' to true tends to
%           reduce the absolute value of correlations.  Trials with too
%           many NaN/outlier values can still be excluded using NaNThresh.
%
%      'NaNThresh',Percent - If the percentage of voxels that are
%           NaNs/outliers exceeds Percent, SimilarityVec will contain NaNs.  
%           This functionality is intended to serve as a safety check to 
%           make sure that results are not driven by removing too many voxels
%           when replaceNaN is (false/default) or driven by replacing too many
%           voxels with mean when replaceNaN is defined as true.  NaNTresh's behavior
%           depends on 'replaceNaN' - when it is true only TRIAL PAIRS with too
%           many NaNs will be considered NaN, when replaceNaN is false,  
%           SimilarityVec will contain only NaNs.
%
% USEAGE:
%[SimilarityMatrix] = MPSMAT(TrialByVoxel1,TrialByVoxel2,'outlier',num,'replaceNaN',bool)

% Check inputs
if ~isa(TrialByVoxel1,'numeric')||~isa(TrialByVoxel2,'numeric')
    error('TrialByVoxel must be a numeric')
end

if ~isempty(varargin)&&all(~strcmpi(varargin,'outlier'))&&all(~strcmpi(varargin,'replaceNaN'))&&all(~strcmpi(varargin,'NaNThresh'))
    error('Unrecognized name-value pair')
end
if ~isempty(varargin)&&any(strcmpi(varargin,'outlier'))
    if isa(varargin{find(strcmpi(varargin,'outlier'))+1},'numeric')
        out = varargin{find(strcmpi(varargin,'outlier'))+1};
    else
        error('Outlier option must be followed by a numeric')
    end
else
    out = inf;
end
if ~isempty(varargin)&&any(strcmpi(varargin,'replaceNaN'))
    if logical(varargin{find(strcmpi(varargin,'replaceNaN'))+1})||any(varargin{find(strcmpi(varargin,'replaceNaN'))+1}==[0 1])
        repnan = varargin{find(strcmpi(varargin,'replaceNaN'))+1};
    else
        error('replaceNaN option must be followed by a logical or 1/0')
    end
else
    repnan = false;
end
if ~isempty(varargin)&&any(strcmpi(varargin,'NaNThresh'))
    if isa(varargin{find(strcmpi(varargin,'NaNThresh'))+1},'numeric')&&varargin{find(strcmpi(varargin,'NaNThresh'))+1}>=0&&varargin{find(strcmpi(varargin,'NaNThresh'))+1}<=100
        percent = varargin{find(strcmpi(varargin,'NaNThresh'))+1};
    else
        error('NaNThresh option must be followed by a numeric between 0 and 100')
    end
else
    percent = 100;
end

dim_tbv1 = size(TrialByVoxel1);
dim_tbv2 = size(TrialByVoxel2);
if dim_tbv1(2)~=dim_tbv2(2)
    error('Data blocks must have same number of voxels')
end
if dim_tbv1(1)~=dim_tbv2(1)
    error('Data blocks must have same number of trials')
end
%% Main function section

%Remove Outliers (these ugly one liners remove outliers calculated across trials)
TrialByVoxel1(bsxfun(@ge,abs(bsxfun(@minus,TrialByVoxel1,nanmean(TrialByVoxel1,2))),out*nanstd(TrialByVoxel1)))=NaN;
TrialByVoxel2(bsxfun(@ge,abs(bsxfun(@minus,TrialByVoxel2,nanmean(TrialByVoxel2,2))),out*nanstd(TrialByVoxel2)))=NaN;
% Replace NaN
if repnan
    [I,J] = ind2sub(size(TrialByVoxel1),find(isnan(TrialByVoxel1)));
    logical_tbv1 = zeros(dim_tbv1);%This logical array is to make sure NaNTresh will properly remove trials
    for i = 1:length(I)
        TrialByVoxel1(I(i),J(i)) = nanmean(TrialByVoxel1(I(i),:),2);
        logical_tbv1(I(i),J(i))=1;
    end
    
    [I,J] = ind2sub(size(TrialByVoxel2),find(isnan(TrialByVoxel2)));
    logical_tbv2 = zeros(dim_tbv2);
    for i = 1:length(I)
        TrialByVoxel2(I(i),J(i)) = nanmean(TrialByVoxel2(I(i),:),2);
        logical_tbv2(I(i),J(i))=1;
    end
    % check percentage parameter if fails single trials are removed
    TrialByVoxel1(sum(logical_tbv1,2)./dim_tbv1(2)>=(percent/100),:)=NaN;
    TrialByVoxel2(sum(logical_tbv2,2)./dim_tbv2(2)>=(percent/100),:)=NaN;
else % remove voxels
    [I,J] = ind2sub(size(TrialByVoxel1),find(isnan(TrialByVoxel1)|isnan(TrialByVoxel2)));
    logical_tbv1 = zeros(dim_tbv1);
    logical_tbv2 = zeros(dim_tbv2);
    
    for i = 1:length(I)
        logical_tbv1(I(i),J(i))=1;
        logical_tbv2(I(i),J(i))=1;
    end
    J = unique(J);
    for i = 1:length(J)
        TrialByVoxel1(:,J(i)) = [];
        TrialByVoxel2(:,J(i)) = [];
        J = J-1;
    end
    % check percentage parameter if fails all trials removed
    if any(sum(logical_tbv1,2)./dim_tbv1(2)>=(percent/100))||any(sum(logical_tbv2,2)./dim_tbv2(2)>=(percent/100))
        SimilarityVec = NaN(dim_tbv1(1),1);
        return;
    end
end

% Perform correlations
TrialByVoxel1 = bsxfun(@minus,TrialByVoxel1,nanmean(TrialByVoxel1,2));
TrialByVoxel2 = bsxfun(@minus,TrialByVoxel2,nanmean(TrialByVoxel2,2));

SimilarityVec=zeros(size(TrialByVoxel1,1),1);
for i = 1:size(TrialByVoxel1,1)
    SimilarityVec(i,1) = (TrialByVoxel1(i,:)*TrialByVoxel2(i,:)')./(sqrt(TrialByVoxel1(i,:)*TrialByVoxel1(i,:)').*sqrt(TrialByVoxel2(i,:)*TrialByVoxel2(i,:)'));
end

end

