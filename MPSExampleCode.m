%% An example multivariate pattern similarity analysis script.
% (also known as representational similarity) 
% This code demonstrates the process of performing MPS analysis, similar to:
%
% Complementary roles of human hippocampal subregions during retrieval of 
%   spatiotemporal context. (2014) MS Copara, AS Hassan, CT Kyle, LA Libby, 
%   C Ranganath, AD Ekstrom. Journal of Neuroscience 34 (20), 6834-6842
%
% Stokes, J., Kyle, C., & Ekstrom, A. D. (2015). Complementary Roles of 
%   Human Hippocampal Subfields in Differentiation and Integration of 
%   Spatial Context. Journal of Cognitive Neuroscience.
%
% Successful retrieval of competing spatial environments in humans involves 
%   hippocampal pattern separation mechanisms. (2015). CT Kyle, JD Stokes, 
%   JS Lieberman, AS Hassan, AD Ekstrom Elife 4, e10499
%
% Kyle, C., Smuda, D., Hassan, A., & Ekstrom, A. (2015). Roles of human 
%   hippocampal subfields in retrieval of spatial and temporal context. 
%   Behavioral Brain Research
%
% --- REQUIRES SPM IN PATH ---
% --- REQUIRES THAT ADVANCED NORMALIZATION TOOLS IS INSTALLED ---
%% General Environment Variables
% set main directory 
main_dir = ['/Users/colin/Documents/MATLAB/Imaging/MPS_Example/'];
% Masks are 3D binary nifti images in the same space as your Betas
mask_dir = [main_dir,'Masks/'];
% Betas are 3D .hdr/.img images (output from an SPM single trial model)
beta_dir = [main_dir,'BetasMultiHRF/'];
% This is where we'll save the Region of Interest analysis results
roi_dir = [main_dir,'ROIData/'];
% This is where we'll save the Searchlight analysis results
search_dir = [main_dir,'SearchData/'];

% set subject names etc
% data is stored in individual subject folders:
% i.e., mask_dir/Subject01/masks.nii
%       beta_dir/Subject05/beta_0001.hdr
subnums = {'02'; '03'; '04'; '05'; '06'; '07'; '08'; '09';...
    '10'; '11'; '12'; '13'; '14'; '15'; '16'; '17'; '18'; '19'; '20';};
for i = 1:length(subnums)
    subjects{i,1} = ['Subject',subnums{i},'/'];
    subjectsn{i,1} = ['Subject',subnums{i}];
end

% mask strings follow naming convention:
%   mask_dir/Subject01/coregdSubject01_Left_CA1.nii
%   mask_dir/Subject01/coregdSubject01_Right_CA1.nii
%   mask_dir/Subject05/coregdSubject05_Left_CA1.nii
%   mask_dir/Subject05/coregdSubject05_Right_CA1.nii
masks = {'_Left_CA1.nii'; '_Left_CA3.nii'; ...
    '_Left_CAant.nii'; '_Left_SUB.nii'; '_Left_PHC.nii';...
    '_Left_EC.nii'; '_Left_PC.nii';...
    '_Right_CA1.nii'; '_Right_CA3.nii'; ...
    '_Right_CAant.nii'; '_Right_SUB.nii'; '_Right_PHC.nii';...
    '_Right_EC.nii'; '_Right_PC.nii';};
% for plotting
maskstrs = {'LCA1'; 'LCA3'; 'LCAant'; 'LSUB'; 'LPHC'; 'LEC'; 'LPC';...
    'RCA1'; 'RCA3'; 'RCAant'; 'RSUB'; 'RPHC'; 'REC'; 'RPC';};

% We're going to assume the data of Kyle et al., 2015 in eLife
% In this task subjects retrieved memories from 1 of 4 cities, two times
% each.  Therefore, there are 8 retrieval/recording blocks.  Subjects were asked
% different questions while retrieving a city for the second time so I coded 
% if up as if there are 8 cities.  However city5 is the repeat of city1,
% city 6 the repeat of city 2, etc.
% This also means I have one city condition per retrieval
% block so this variable will also serve as a "recording block" variable.  Keeping
% different retrieval blocks seperate is important for MPS analysis because
% correlations within block will be artificially high and correlations
% between blocks will be lower.  No analysis should never mix within and 
% between blocks in the same analysis/test

% condition/city strings, 8 retrieval blocks but only 4 cities/conditions,
% they're named cities 1-8 out of convenience.
% city 5 is a repeat of city 1, city 6 is a repeat of city 2, etc...
for i = 1:8
    cities{i,1} = ['city',num2str(i),'/'];
    city{i,1} = ['city',num2str(i)];
end
% extension strings for the lazy
nii = '.nii';
mat = '.mat';
% functional image size
func_dim = [118 120 36];%Put the 3D dimensions of your functional images here


%% ROI MPS Analysis (Region of Interests defined by masks 
% (mine are hippocampal subfields CA1, CA3, ...)
%
% Final data product is MPS measurements, with 1 value per
% subject/region/condition
ROI_MPS =1;%toggle analysis on/off
if ROI_MPS
    for sub = 1:length(subjects)
        % first lets load the masks
        cd([mask_dir,subjects{sub}])
        % generate mask names
        for mas = 1:numel(masks)
            masknames{mas} = strcat('coregd',subjectsn{sub},masks{mas});
        end
        %load them
        [maskdata.sub(sub).maskinds,~,maskdata.sub(sub).hdr] ...
            = getmasks(masknames);
        %% LOAD BETAS FROM EACH ROI
        % next lets load data, in this analysis, each block of data is from a
        % different city (ie condition), this section may look a little
        % different if you have multiple conditions per block.
        
        for block = 1:numel(cities)
            cd([beta_dir,subjects{sub},cities{block}])
            %% MAKE CONDITION/BETA NAME VARIABLE
            % Things to consider... does order matter? Do I eventually want to
            % make correlations between specific trials? Or do I just want to
            % make correlations between all pairs within or between conditions?
            %  Do I want to just re-order later?
            
            % For now lets not pay attention to order, but I might do the
            % following if I do want to do a specific order here:
            % ordercond1 = [18, 5, 1, 7, ...];
            % ordercond2 = [20, 16, 6, 3, ...];
            % for i = 1:numel(ordercond1)
            %   betacondition1{i} = ['beta_',num2str(ordercond1(i),'%04i'),'.hdr'];
            % end ...
            
            % Just grab all betas with no regard to order
            for i = 1:20
                betascondition1{i} = ['beta_',num2str(i,'%04i'),'.hdr'];
            end
            
            % I only have one condition per block but you'll want to
            % create different beta condition strings for each condition here
            %Example
            %       for i = 1:numel(ordercond1)
            %           betascondition2{i} = ['beta_',num2str(ordercond2(i),'%04i'),'.hdr'];
            %       end
            
            % load betas values contained within each mask for all trials of a
            % condition:
            ROIdata.sub(sub).block(block).condition1 = getdata3d(betascondition1,maskdata.sub(sub).maskinds);
            % data is formatted TrialByVoxel{mas}
        end
        %% Normally I'd just finish this loop and save all the data as a .mat
        % to make for easier loading later, but lets go ahead and do our MPS
        % analysis right here for the sake of this tutorial
        
        % A COMPLETE MPS ANALYSIS HAPPENS ONCE PER MASK!
        for mas = 1:numel(masks)
            % I've made two basic MPS functions.  They both have similar inputs but
            % make different types of correlations.  The first is mpsmat:
            % [SimilarityMatrix] = mpsmat(TrialByVoxel1,TrialByVoxel2)
            % This function makes correlations between all pairwise combinations of
            % trials.
            % The second is mpsvec:
            % [SimilarityVec] = mpsvec(TrialByVoxel1,TrialByVoxel2)
            % This function makes correlations between matched pairs of Trials.
            % ex: TrialByVoxel1(1,:) vs. TrialByVoxel2(1,:), TrialByVoxel1(2,:) vs.
            % TrialByVoxel2(2,:), etc ...
            
            [SimilarityMatrix] = mpsmat(ROIdata.sub(sub).block(1).condition1{mas},ROIdata.sub(sub).block(5).condition1{mas});
            [SimilarityVec] = mpsvec(ROIdata.sub(sub).block(1).condition1{mas},ROIdata.sub(sub).block(5).condition1{mas});
            
            if(~exist('MakePlots','var'))
                MakePlots = input('Do you want to plot the data? (1-yes/0-no):  ');
            end
            
            if MakePlots
                figure
                imagesc(SimilarityMatrix);
                
                title(['MPS Matrix ',maskstrs{mas},' ',subjectsn{sub}])
                figure
                imagesc(SimilarityVec);
                title(['MPS Vector ',maskstrs{mas},' ',subjectsn{sub}])
                
                input('Press any Button to Continue with the next mask ...')
            end
            
            % Note that SimilarityVec is the diagonal of SimilarityMatrix
            % One MPS value per subject, mask, condition
            % where conditions are fullmatrix MPS and vec MPS in this case
            MPSDATA.mask(mas).fullmatrix(sub,1)=nanmean(nanmean(SimilarityMatrix));
            MPSDATA.mask(mas).vec(sub,1)=nanmean(SimilarityVec);
        end
        
        
        
        
        
    end
    cd(mask_dir)
    save('mask_inds.mat','-struct','maskdata')
    mkcd(roi_dir)
    save('roi_data_raw.mat','-struct','ROIdata')
    save('roi_data_group.mat','-struct','MPSDATA')
end

%% Searchlight MPS Analysis
% This analysis is more complicated.
% 1. Run searchlight (run analysis in convolution through brain) in each
%   subject. Output is Statistical image with MPS value at each location.
% 2. Warp statistical images into common space across subjects.  This uses
%   ANTS to perform non-rigid warping to template subject.
% 3. Run group analysis 
%   A. Run contrast test at each voxel.
%   B. Cluster results to correct for multiple comparisons.
%
% Result is significant clusters where hypothesis is true.
SEARCH_MPS = 1;
if SEARCH_MPS
    % load masks Data
    cd(mask_dir)
    maskinds = load('mask_inds.mat');
    
    for sub = 1:length(subjects)% this loop takes ~20 seconds per loop, will
        % scale up as MPS conditions and comparisons are added
        
        %% create search space!!
        % the process is similar to getting the masks in the ROI_MPS
        % section, however now we're getting one mask per searchlight in
        % the entire MTL
        
        % first we create a searchlight shape
        n=linspace(-1,1,5);
        [x,y,z]=ndgrid(n,n,n);
        search_sphere=(sqrt(x.^2 + y.^2 + z.^2)<=1);
        % just made a sphere, lets snip the tips
        search_sphere(:,:,1) = [];
        search_sphere(:,:,end) = [];
        % Ok next function does a lot of magic
        % first, it takes the mask_inds you give it and creates a solid
        % volume that encompasses them, then it moves the search_sphere
        % throughout the volume to check all the places it can
        % occupy, then it returns those places it can occupy in a cell that
        % can be used synonymously to the output of getmasks
        % One slight difference is that we need to keep the "searchlight
        % centers" in order to know where to place our statistics for each
        % searchlight volume later
        [searchmaskdata.sub(sub).maskinds,searchmaskdata.sub(sub).centers] = ...
            getsearchinds(func_dim,maskinds.sub(sub).maskinds,search_sphere,2);
        
        
        % Next we need to collect the parameter estimates from each
        % searchcenter just like we did with the masks in ROI_MPS
        for block = 1:numel(cities)
            cd([beta_dir,subjects{sub},cities{block}])
            %% MAKE CONDITION/BETA NAME VARIABLE
            % Things to consider... does order matter? Do I eventually want to
            % make correlations between specific trials? Or do I just want to
            % make correlations between all pairs within or between conditions?
            %  Do I want to just re-order later?
            
            % For now lets not pay attention to order, but I might do the
            % following if I do want to do a specific order here:
            % ordercond1 = [18, 5, 1, 7, ...];
            % ordercond2 = [20, 16, 6, 3, ...];
            % for i = 1:numel(ordercond1)
            %   betacondition1{i} = ['beta_',num2str(ordercond1(i),'%04i'),'.hdr'];
            % end ...
            
            % Just grab all betas with no regard to order
            for i = 1:20
                betascondition1{i} = ['beta_',num2str(i,'%04i'),'.hdr'];
            end
            
            % I only have one condition per block but you'll want to
            % create different beta condition strings for each condition here
            %Example
            %       for i = 1:numel(ordercond1)
            %           betascondition2{i} = ['beta_',num2str(ordercond2(i),'%04i'),'.hdr'];
            %       end
            
            % load betas values contained within each mask for all trials of a
            % condition:
            searchdata.sub(sub).block(block).condition1 = getdata3d(betascondition1,searchmaskdata.sub(sub).maskinds);
            %data is formatted TrialByVoxel{mas}
        end
        
        % Here's where the searchlight happens, we loop through each
        % searchlight location and perform the same analysis we did in the
        % ROI analysis above
        % Create image that will serve as statistical map
        statimg1 = NaN(func_dim);
        statimg2 = NaN(func_dim);
        
        for search = 1:numel(searchmaskdata.sub(sub).centers)
            % I've made two basic MPS functions.  They both have similar inputs but
            % make different types of correlations.  The first is mpsmat:
            % [SimilarityMatrix] = mpsmat(TrialByVoxel1,TrialByVoxel2)
            % This function makes correlations between all pairwise combinations of
            % trials.
            % The second is mpsvec:
            % [SimilarityVec] = mpsvec(TrialByVoxel1,TrialByVoxel2)
            % This function makes correlations between matched pairs of Trials.
            % ex: TrialByVoxel1(1,:) vs. TrialByVoxel2(1,:), TrialByVoxel1(2,:) vs.
            % TrialByVoxel2(2,:), etc ...
            
            [SimilarityMatrix] = mpsmat(searchdata.sub(sub).block(1).condition1{search},...
                searchdata.sub(sub).block(5).condition1{search},'NaNThresh',70);
            [SimilarityVec] = mpsvec(searchdata.sub(sub).block(1).condition1{search},...
                searchdata.sub(sub).block(5).condition1{search},'NaNThresh',70);
            
            
            % Note that SimilarityVec is the diagonal of SimilarityMatrix
            
            % Save in statistical map for warping and group analysis
            statimg1(searchmaskdata.sub(sub).centers(search)) = nanmean(nanmean(SimilarityMatrix));
            statimg2(searchmaskdata.sub(sub).centers(search)) = nanmean(SimilarityVec);
            
        end
        %% Multiply stat images by 10000, ANTS uses single precision so small
        % Decimals will become zero!
        statimg1=10000*statimg1;
        statimg2=10000*statimg2;
        
        % Get hdr that has the same affine matrix as masks and betas
        cd([mask_dir,subjects{sub}])
        hdr = spm_vol(['coregd',subjectsn{sub},nii]);
        % create and move to storage directory
        mkcd([search_dir,subjects{sub}])
        %specify name of statistical map image
        hdr.fname = 'statimg1.nii';
        %write statistical image
        spm_write_vol(hdr,statimg1);
        %specify name of statistical map image
        hdr.fname = 'statimg2.nii';
        %write statistical image
        spm_write_vol(hdr,statimg2);
        
    end
end

%% Warp to template subject
WARP = 1;
if WARP
    % load mask data
    cd(mask_dir)
    maskinds = load('mask_inds.mat');
    GENWARPMASKS = 1;
    if GENWARPMASKS
        %% Create the guide mask used for image registration
        % image registration weights guide images created from the masks, 
        % and the structural images in order to produce best results.
        % Guide masks have intensity values highest in the centroid of the
        % mask which fall away to the ends of the mask.  This intensity
        % gradient ensures smooth warping and maintains topology
        
        maskweights = [500,100,300,200,50,80,30,500,100,300,200,50,80,30];
        for sub = 1:length(subjects)
            cd([mask_dir,subjects{sub}])
            normalmask('Guide.nii',maskinds.sub(sub).maskinds,maskweights,maskinds.sub(sub).hdr);
        end
    end
    
    %% Generate warping parameters and do warp (once you have the warping 
    % pamaters estimated once, there's no need to do it again,
    % comment out the estimate warp line if you want to go back and )
    DOWARP = 1;
    if DOWARP
        for sub=1:length(subjects)
            %% CD to mask dir so that full paths are only necessary for the
            %template subject
            %Estimate warping parameters using estimatewarp (and ANTS)
            cd([mask_dir,subjects{sub}])
            estimatewarp('estwarp','MSQ',[mask_dir,subjects{5},'Guide.nii'],'Guide.nii',.6,0,...
                'CC',[mask_dir,subjects{5},'coregd',subjectsn{5},nii],['coregd',subjectsn{sub},nii],.4,5);
            %% Now apply warp to all masks/structurals
            submasks = file('coregd*');
            mkdir('warpedtotemplate');
            for mas = 1:length(masks)
                applywarp('estwarp',submasks{mas},['warpedtotemplate/',submasks{mas}],[mask_dir,subjects{5},'Guide.nii'])
            end
            %% Now apply warp to all searchlight images
            cd([search_dir,subjects{sub}])
            statimgs = file('stat*');
            mkdir('warpedtotemplate/')
            for st = 1:length(statimgs)
                applywarp([mask_dir,subjects{sub},'estwarp'],statimgs{st},['warpedtotemplate/',statimgs{st}],[mask_dir,subjects{5},'Guide.nii']);
            end
        end
    end
end

%% Group Analysis
GROUP = 0;
if GROUP
    %% Load each subject's groupspace Data!
    for sub=1:length(subjects)
        cd([search_dir,subjects{sub},'warpedtotemplate/'])
        statimgs = file('stat*');
        %%first get valid inds for each subject
        data = cell2mat(getdata3d(statimgs));
        validinds{sub} = find(data(1,:)~=0);%ANTS turns NaNs to zeros in the warping process
        for i = 2:length(statimgs)
            validinds{sub} = intersect(validinds{sub},find(data(i,:)~=0));
        end
    end
    %% Now we just want the inds common to all subjects
    groupinds = validinds{1};
    for sub=2:length(subjects)
        groupinds = intersect(groupinds,validinds{sub});
    end
    % Now group inds is one big mask, repeat process but just the inds we
    % want
    for sub = 1:length(subjects)
        cd([search_dir,subjects{sub},'warpedtotemplate/'])
        statimgs = file('stat*');
        for i = 1:length(statimgs)
            GroupData.(statimgs{i}(1:(strfind(statimgs{i},'.')-1)))(sub,:) =...
                cell2mat(getdata3d(statimgs{i},groupinds));
        end
    end
    
    %% Now we can do our contrasts
    [h,p,c,t] = ttest(GroupData.statimg1,GroupData.statimg2);
    contrastIMG = zeros(func_dim);
    contrastIMG(groupinds)=t.tstat;
    %% Now we cluster!
    [Clusters] = clusterimage(contrastIMG,tinv(.95,length(subjects)));
    
    %% Plot largest cluster
    % Get template subject's structural
    hdr = spm_vol([mask_dir,subjects{5},'coregd',subjectsn{5},nii]);
    templateBrain = spm_read_vols(hdr);
    templateBrain(Clusters{1})=max(templateBrain(:));
    figure
    for i=1:length(Clusters{1})
        clinds(i) = find(groupinds==Clusters{1}(i));
    end
    Data = [nanmean(GroupData.statimg1(:,clinds),2),nanmean(GroupData.statimg2(:,clinds),2)];
    subplot(2,1,2)
    errorplot(Data)
    [I,J,K] = ind2sub(func_dim,Clusters{1});
    subplot(2,1,1)
    for i = min(K):max(K)
        imagesc(imrotate(templateBrain(:,:,i),90))
        colormap('gray')
        drawnow
        pause(.5)
    end
    
end