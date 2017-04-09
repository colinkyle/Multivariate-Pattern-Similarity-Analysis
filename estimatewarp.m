function [ output_args ] = estimatewarp(warpprefix,measure,fixed,moving,weight,parameter,varargin)
%ESTIMATEWARP this is really just some code that makes Advanced
% Normalization Tools (ANTS) a little easier to call from MATLAB
% USE:
% estimatewarp(warpprefix,measure,fixed,moving,weight,parameter)
%
% if you get something like "/bin/bash: ANT: command not found" it means
% that you haven't added ANTS to your local profile, to do so run the
% following in matlab:
%   edit ~/.bash_profile
%   
%   #ANTs setup
%   ANTSPATH=/Users/colin/Applications/MyCode/antsbin/bin
%   PATH=${PATH}:${ANTSPATH}
%   Then close matlab and open from terminal using something like:
%   /Applications/MATLAB_R2015a.app/bin/matlab
%
%   Optionally, you can add as many measures with different images as you want:
%   [ output_args ] =
%   estimatewarp(warpout,measure1,fixed1,moving1,weight1,parameter1,...
%        measure2,fixed2,moving2,weight2,parameter2,...);

% Dummy variables

meas{1} = measure;
fix{1} = fixed;
mov{1} = moving;
wei{1} = weight;
param{1} = parameter;

% Check varargin parameters
if mod(numel(varargin),5)
    error('You didn''t use the right number of inputs, see help')
end
ct=2;
for i = 1:5:numel(varargin)
    i
    if isa(varargin{i},'char')
        meas{ct} = varargin{i};
    else
        disp(varargin{i})
        error('measure must be a char')
    end
    
    if isa(varargin{i+1},'char')
        fix{ct} = varargin{i+1};
    else
        disp(varargin{i+1})
        error('fixed must be a char')
    end
    
    if isa(varargin{i+2},'char')
        mov{ct} = varargin{i+2};
    else
        disp(varargin{i+2})
        error('moving must be a char')
    end
    
    if isa(varargin{i+3},'numeric')&&varargin{i+3}<=1
        wei{ct} = varargin{i+3};
    else
        disp(varargin{i+3})
        error('weight must be a numeric less than or equal to 1')
    end
    
    if isa(varargin{i+4},'numeric')
        param{ct} = varargin{i+4};
    else
        disp(varargin{i+3})
        error('parameter must be a numeric less than or equal to 1')
    end
    ct=ct+1;
end

callstr = '!ANTS 3 ';
for i = 1:numel(meas)
    callstr = [callstr,'-m ',meas{i},'['];
    callstr = [callstr,fix{i},', '];
    callstr = [callstr,mov{i},','];
    callstr = [callstr,num2str(wei{i}),','];
    callstr = [callstr,num2str(param{i}),'] '];
end
%--number-of-affine-iterations 5000x5000x5000 --affine-gradient-descent-option 0.5x0.9x1.e-4x1.e-4
callstr = [callstr,'--number-of-affine-iterations 0 -i 5x5x1 -o '];
%callstr = [callstr,' -o '];
callstr = [callstr,warpprefix,'.nii'];
callstr = [callstr,' -t SyN[.5] -r Gauss[1.5,1.5] ']
%callstr = [callstr,' -t Elast[1] -r Gauss[0.5,3] ']

eval(callstr)

end

