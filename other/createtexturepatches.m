%% init vars
bverbose = true;
normmethod = 'zscore'; % normalization method: 'zscore','wm'
testAUC = [];
allData = {};
volumes = {};
masks = {};
allGroundTruth = {};
patients = [1]
% patients = [1,2,4,5,6,7,8,10,12,14,15,16,17,18,19]; %no AUC = 0 no disease progression from combine.m
for J = 1:length(patients)
    num = num2str(patients(J));
    directoryName = strcat('../genMSdata2/patient', num);
    if ~exist(directoryName, 'dir')
        mkdir(directoryName)
    end
end
%% load data from nii into allData and allGroundTruth
for J = 1:length(patients)
    patient_number = patients(J); %change the particular patient's number
    num = num2str(patient_number);
    folder = strcat('../MSpatientdata/patient', num);
    % Get a list of all files in the folder with the desired file name pattern
    filePattern = fullfile(folder, '*.nii*'); % Change to whatever pattern you need.
    theFiles = dir(filePattern);
    s.patient_number = num;

    for k = 1 : length(theFiles)
        baseFile = theFiles(k).name;
        fullFileName = fullfile(folder, baseFile);
        fprintf(1, 'Now reading %s\n', fullFileName);
        %     variable allocation
        if strfind(baseFile, '1_T1')
            s.t1_s1file = fullFileName;
        elseif strfind(baseFile, '1_T2')
            s.t2_s1file = fullFileName;
        elseif strfind(baseFile, '1_FLAIR')
            s.flair_s1file = fullFileName;
        elseif strfind(baseFile, '2_T1')
            s.t1_s2file = fullFileName;
        elseif strfind(baseFile, '2_T2')
            s.t2_s2file = fullFileName;
        elseif strfind(baseFile, '2_FLAIR')
            s.flair_s2file = fullFileName;
        elseif strfind(baseFile, 'gt3')
            s.gtfile = fullFileName;
        elseif strfind(baseFile, 'mask')
            s.maskfile = fullFileName;
        end
    end
    %% load and normalize data

    nii = load_nii(s.maskfile);
    mask = nii.img;
    mask = logical(mask);
    nii = load_nii(s.gtfile);
    % load the new gt:
    gt = getfield(load_nii([folder,'/patient', num,'_gt3.nii']),'img');
    fields = fieldnames(s);
    fields = setdiff(fields,{'maskfile','gtfile','patient_number'});%%skip patient number, also skip the brain mask and the ground truth mask
    ints3d = struct;

    if strcmp(normmethod,'zscore')
        for N = 1:numel(fields)
            field = fields{N};
            value = getfield(s, field);
            nii = load_nii(value);
            im = nii.img;
            r = strrep(field,'file','');
            temp = im(mask);
            centered = (temp - mean(temp)) ;
            %         sum = 0;
            %         for i = 1:numel(centered)
            %             sum = sum + (centered(i, 1)^2);
            %         end
            stddev = std(double(temp));
            %     mean(temp)
            %     ints.(strcat(r,'ints')) = centered; %0 centered, now all stored in ints - should this be before or after difference calculation?
            ints3d.(strcat(r,'ints')) = zeros(size(mask));
            ints3d.(strcat(r,'ints'))(mask) = centered/stddev;
        end
       %% normalize the intensities wrt WM
    elseif strcmp(normmethod,'wm')
        [ints3d.t1_s1ints, ints3d.t2_s1ints, ints3d.flair_s1ints] = normalizeWM(getfield(load_nii(s.t1_s1file),'img'),...
            getfield(load_nii(s.t2_s1file),'img'),...
            getfield(load_nii(s.flair_s1file),'img'), mask);
    end
   
    %% calc distances
    dist.t1_dist = ints3d.t1_s2ints - ints3d.t1_s1ints;
    dist.flair_dist = ints3d.flair_s2ints - ints3d.flair_s1ints;
    dist.t2_dist = ints3d.t2_s2ints - ints3d.t2_s1ints;
    %% only hyperintense t2 voxel selection mask
%     vsmask = zeros(size(ints3d.t2_s1ints));
%     vsmask(dist.flair_dist>std(dist.flair_dist(mask))) =1;
%     vsmask = logical(vsmask);
    
    temp = zeros(size(dist.flair_dist));
    vsmask = zeros(size(dist.flair_dist));
    for i=1:size(dist.flair_dist,3)
        temp(:,:,i) = medfilt2(dist.flair_dist(:,:,i));
    end
    vsmask(temp>(std(temp(:)))) = 1;
    vsmask = logical(vsmask);
    notzero = nnz(vsmask)
    %% additional info to save
    labels = gt(vsmask);
    inds = find(vsmask==1);

    %% get indices for each selected voxel
%     inds = reshape(IC,size(vsmask,1),size(vsmask,2),size(vsmask,3)); 
%     inds(vsmask==0) = [];
%     inds = inds';
%     size(inds)
    %% balanced voxel selection mask based on t2 intensities
%     voxel_selection_mask = zeros(size(gt));
%     voxel_selection_mask(gt==1) = 1; % select all positive class observations
%     voxel_selection_mask(randsample(find(mask>0 & gt~=1 & ints3d.t2_s1ints>0),nnz(gt==1))) = 1; % randomly select the same number of other observations
%     voxel_selection_mask = logical(voxel_selection_mask);
        %% calculate texture features
    textureProps ={};
    
    IMG = {ints3d.t1_s1ints,ints3d.t2_s1ints,ints3d.flair_s1ints,...
        dist.t1_dist, dist.t2_dist, dist.flair_dist,...
        ints3d.t1_s2ints,ints3d.t2_s2ints,ints3d.flair_s2ints};
    for ii = 1:numel(IMG)% over all "channels"
        %     vsmtextureProps =[];
        tmpTextureProps =[];
        patchSize = 21;
        empty = ones(patchSize);
        %     empty((patchSize * patchSize + 1)/2) = 1;
        for i=1:size(IMG{ii},3) %by slice
            
            patches = nlfilter(IMG{ii}(:,:,i), [patchSize patchSize], @(block) {block});
            
            %         vsmpatches = patches;
            
            %on mask
            %             maskslice = mask(:,:,i);
            %             patches(maskslice == 0) = [];
            %             for k=1:numel(patches)
            %                 temp = patches{k};
            % %                 temp = reshape(patches(:,k), patchSize, patchSize);
            %                 [feat,featureGroupID, featureGroupList, featureList] = computeFeats(temp, empty==1);
            %                 textureProps = [textureProps, feat];
            %             end
            
            
            %on voxel selection mask
            %             maskslice = voxel_selection_mask(:,:,i);
            %             vsmpatches(maskslice == 0) = [];
            
            %             for k=1:numel(vsmpatches)
            %                 temp = vsmpatches{k};
            % %                 temp = reshape(patches(:,k), patchSize, patchSize);
            %                 [feat,featureGroupID, featureGroupList, featureList] = computeFeats(temp, empty==1);
            %                 vsmtextureProps = [vsmtextureProps, feat];
            %             end
            maskslice = vsmask(:,:,i);
            patches(maskslice == 0) = [];
            parfor k=1:numel(patches)
                temp = patches{k};
                %                 temp = reshape(patches(:,k), patchSize, patchSize);
                [feat,featureGroupID, featureGroupList, featureList] = computeFeats(temp, empty==1);
                tmpTextureProps = [tmpTextureProps, feat];
            end
            
            %             %on voxel selection mask
            %             maskslice = voxel_selection_mask(:,:,i);
            %             vsmpatches(maskslice == 0) = [];
            %             for k=1:numel(vsmpatches)
            %                 temp = vsmpatches{k};
            % %                 temp = reshape(patches(:,k), patchSize, patchSize);
            %                 [feat,featureGroupID, featureGroupList, featureList] = computeFeats(temp, empty==1);
            %                 vsmtextureProps = [vsmtextureProps, feat];
            %             end
        end
        textureProps{ii} = tmpTextureProps;
        size(textureProps{ii})
        textureProps{ii} = textureProps{ii}';
        %     vsmtextureProps = vsmtextureProps';
        %     testsize = size(vsmtextureProps)
    end
     %%    save textureprops

    directoryName = strcat('../MSproject/texturepatches9ch/patient', num);
    if ~exist(directoryName, 'dir')
        mkdir(directoryName)
    end

%     h5create(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'), '/dataset1', size(textureProps));
%     h5write(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'), '/dataset1', textureProps)
%     h5create(strcat(directoryName, '/patient', num,'_balancedTextureProps.h5'), '/dataset1', size(vsmtextureProps));
%     h5write(strcat(directoryName, '/patient', num,'_balancedTextureProps.h5'), '/dataset1', vsmtextureProps)

    channelnames = {'/t1_s1','/t2_s1','/flair_s1',...
        '/delta_t1','/delta_t2','/delta_flair',...
        '/t1_s2','/t2_s2','/flair_s2'};
    for ii = 1:numel(channelnames)
        h5create(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'),...
            channelnames{ii}, size(textureProps{ii}));
        h5write(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'), ...
            channelnames{ii}, textureProps{ii})
    end
    
    h5create(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'),...
        '/labels',size(labels));
    h5write(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'), ...
        '/labels',labels)
    h5create(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'),...
        '/inds',size(inds));
    h5write(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'), ...
        '/inds', inds)
%     h5create(strcat(directoryName, '/patient', num,'_balancedTextureProps.h5'), '/dataset1', size(vsmtextureProps));
%     h5write(strcat(directoryName, '/patient', num,'_balancedTextureProps.h5'), '/dataset1', vsmtextureProps)

end