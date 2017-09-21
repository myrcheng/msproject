% adding texture features/neighborhood analysis
% dependencies: NifTI, GLRL(add to selected folders and subfolders), computeFeats.m
for iteration = 1
%% init vars
bverbose = true;
normmethod = 'zscore'; % normalization method: 'zscore','wm'
testAUC = [];
alltestData = {};
alltrainData = {};
alltrainGroundTruth = {};
alltestGroundTruth = {};
volumes = {};
masks = {};
texturecoef = {};
set1 = [];
set2 = []
seflair = []

semedt1 = [];
semedt2 = []
semedflair = []


% patients = [1]
patients = [1,2,4,5,6,7,8,10,12,14,15,16,17,18,19]; %no AUC = 0 no disease progression from combine.m
%% load data from nii into allData and allGroundTruth
% patients =[1];
for J = 1:length(patients)
    patient_number = patients(J); %change the particular patient's number
    num = num2str(patient_number);
    fprintf(strcat('\n now reading patient', num, '\n'))
    folder = strcat('../MSpatientdata/patient', num);
    % Get a list of all files in the folder with the desired file name pattern
    filePattern = fullfile(folder, '*.nii*'); % Change to whatever pattern you need.
    theFiles = dir(filePattern);
    s.patient_number = num;

    for k = 1 : length(theFiles)
        baseFile = theFiles(k).name;
        fullFileName = fullfile(folder, baseFile);
%         fprintf(1, 'Now reading %s\n', fullFileName);
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
    bmask = nii.img;
    bmask = logical(bmask);
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
            temp = im(bmask);
            centered = (temp - mean(temp)) ;
            %         sum = 0;
            %         for i = 1:numel(centered)
            %             sum = sum + (centered(i, 1)^2);
            %         end
            stddev = std(double(temp));
            %     mean(temp)
            %     ints.(strcat(r,'ints')) = centered; %0 centered, now all stored in ints - should this be before or after difference calculation?
            ints3d.(strcat(r,'ints')) = zeros(size(bmask));
            ints3d.(strcat(r,'ints'))(bmask) = centered/stddev;
        end
       %% normalize the intensities wrt WM
    elseif strcmp(normmethod,'wm')
        [ints3d.t1_s1ints, ints3d.t2_s1ints, ints3d.flair_s1ints] = normalizeWM(getfield(load_nii(s.t1_s1file),'img'),...
            getfield(load_nii(s.t2_s1file),'img'),...
            getfield(load_nii(s.flair_s1file),'img'), bmask);
    end
   
    %% calc distances
    dist.t1_dist = ints3d.t1_s2ints - ints3d.t1_s1ints;
    dist.flair_dist = ints3d.flair_s2ints - ints3d.flair_s1ints;
    dist.t2_dist = ints3d.t2_s2ints - ints3d.t2_s1ints;
    %%%%%%%%%%START
     %% only hyperintense t2 voxel selection mask 
    vsmask = zeros(size(dist.flair_dist));
    vsmask(ints3d.t1_s1ints>0) = 1;
    vsmask = logical(vsmask);
    Se = nnz(vsmask>0 & gt==1)/nnz(gt==1);
    set1 = [set1, nnz(vsmask)];
    
     vsmask = zeros(size(dist.flair_dist));
    vsmask(ints3d.flair_s1ints>0) = 1;
    vsmask = logical(vsmask);
    Se = nnz(vsmask>0 & gt==1)/nnz(gt==1);
    seflair = [seflair, nnz(vsmask)];
    vsmask = zeros(size(dist.flair_dist));
    vsmask(ints3d.t2_s1ints>0) = 1;
    vsmask = logical(vsmask);
    Se = nnz(vsmask>0 & gt==1)/nnz(gt==1);
    set2 = [set2, nnz(vsmask)];
    
    
    temp = zeros(size(dist.t1_dist));
    vsmask = zeros(size(dist.t1_dist));
    for i=1:size(dist.t1_dist,3)
        temp(:,:,i) = medfilt2(dist.t1_dist(:,:,i));
    end
    vsmask(temp>std(temp(:))) = 1;
    Se = nnz(vsmask>0 & gt==1)/nnz(gt==1);
    semedt1 = [semedt1, nnz(vsmask)];
    
    temp = zeros(size(dist.t2_dist));
    vsmask = zeros(size(dist.t2_dist));
    for i=1:size(dist.t2_dist,3)
        temp(:,:,i) = medfilt2(dist.t2_dist(:,:,i));
    end
    vsmask(temp>std(temp(:))) = 1;
    Se = nnz(vsmask>0 & gt==1)/nnz(gt==1);
    semedt2 = [semedt2, nnz(vsmask)];
    
    temp = zeros(size(dist.flair_dist));
    vsmask = zeros(size(dist.flair_dist));
    for i=1:size(dist.flair_dist,3)
        temp(:,:,i) = medfilt2(dist.flair_dist(:,:,i));
    end
    vsmask(temp>std(temp(:))) = 1;
    Se = nnz(vsmask>0 & gt==1)/nnz(gt==1);
    semedflair = [semedflair, nnz(vsmask)];
    
%     se = [se, Se];
    fprintf('\nSe: %.2f\n',Se)
    fprintf('\nMissed: %d out of %d\n',nnz(vsmask==0 & gt==1), nnz(gt==1))
%         fprintf('\nNnz: %.2f\n',nnz(vsmask>0))
end
    %% load textureprops
    directoryName = strcat('texturepatches3/patient', num);
    hinfo = hdf5info(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'))
    data = h5read(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'),'/dataset1');
    labels = h5read(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'),'/labels');
    inds = h5read(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'),'/inds');
    %% get the labels. 
%     labels = gt(vsmask);
%following: or just create the mask from the indices...
%     if size(data,1) > nnz(vsmask)
%         fprintf('too big')
%         data = data(inds(ismember(inds, find(vsmask==1))));
%     elseif size(data,1) < nnz(vsmask)
       vsmask = zeros(size(vsmask));
       vsmask(inds)= 1;
       vsmask = logical(vsmask);
%     end
    vsmask_balanced = zeros(size(vsmask));
    vsmask_balanced(vsmask>0 & gt==1) = 1; % select all positive class observations
    vsmask_balanced(randsample(find(vsmask>0 & gt == 0),1*nnz(vsmask>0 & gt==1))) = 1; % randomly select the same number of other observations
    vsmask_balanced = logical(vsmask_balanced);

    tmp = zeros(size(vsmask));
    tmp(vsmask) = 1:nnz(vsmask);

    inds_balanced = tmp(vsmask & vsmask_balanced);
    data_norm = data;
    data_norm(isnan(data_norm))=0;
    % data_norm = (data_norm-repmat(min(data_norm,[],1),size(data_norm,1),1))...
    %     ./repmat(max(data_norm,[],1)-min(data_norm,[],1),size(data_norm,1),1);
    data_norm = (data_norm-repmat(mean(data_norm,1),size(data_norm,1),1))...
        ./repmat(std(data_norm,[],1)-min(data_norm,[],1),size(data_norm,1),1);

    data_train = data_norm(inds_balanced,:);
    labels_train = labels(inds_balanced,:);
    %%%%% END

    %% only hyperintense t2 voxel selection mask
%     voxel_selection_mask = zeros(size(ints3d.t2_s1ints));
%     voxel_selection_mask(ints3d.t2_s1ints>0) =1;
%     voxel_selection_mask = logical(voxel_selection_mask);
    %% balanced voxel selection mask based on t2 intensities
%     voxel_selection_mask = zeros(size(gt));
%     voxel_selection_mask(gt==1) = 1; % select all positive class observations
%     voxel_selection_mask(randsample(find(mask>0 & gt~=1 & ints3d.t2_s1ints>0),nnz(gt==1))) = 1; % randomly select the same number of other observations
% %     voxel_selection_mask(randsample(find(mask>0 & gt~=1 & ints3d.t2_s1ints>0),100)) = 1; % randomly select the same number of other observations
%     
% voxel_selection_mask = logical(voxel_selection_mask);
        %% calculate texture features
%     textureProps =[];
%     patchSize = 21;
%     empty = ones(patchSize);
% %     empty((patchSize * patchSize + 1)/2) = 1;
%     for i=1:size(dist.flair_dist,3) %by slice
%             patches = nlfilter(dist.flair_dist(:,:,i), [patchSize patchSize], @(block) {block});
%             maskslice = voxel_selection_mask(:,:,i);
%             patches(maskslice == 0) = []; 
% %             [feat,featureGroupID, featureGroupList, featureList] = cellfun(@computeFeats, patches, 'UniformOutput', false);
% %             textureProps = [textureProps, [feat{:}]];
%             for k = 1:numel(patches) 
%                 temp = patches{k};
% %                 [feat,featureGroupID, featureGroupList, featureList] = cellfun(@computeFeats, patches, 'UniformOutput', false);
% %                 textureProps = [textureProps, [feat{:}]];
% %                 glcm = graycomatrix(temp);
% %                 grayprops = graycoprops(glcm);
% % %                 temp = reshape(patches(:,k), patchSize, patchSize);
% %                 [feat,featureGroupID, featureGroupList, featureList] = computeFeats(temp, logical(empty));
%       
% %             end
% %     end   
% %   
%                 [feat,featureGroupID, featureGroupList, featureList] = computeFeats(temp, empty==1);
%                 textureProps = [textureProps, feat];
%             end
%     end
%     size(textureProps)
    %% fix texture properties: elim NaN and normalize
%     textureProps(isnan(textureProps))=0;
%     textureProps = textureProps';
    %% normalize scaling between 0 and 1
%     denominator = max(textureProps) - min(textureProps);
%     sub = min(textureProps);
%% normalize using z-scores
%     denominator = std(textureProps);
%     size(denominator)
%     sub = mean(textureProps);
%     %% appply normalization
%     denominator = repmat(denominator, [size(textureProps,1),1]);
%     sub = repmat(sub, [size(textureProps,1),1]);
%     result = (textureProps - sub)./denominator;
% 
%     size(result)
% %% find delta_flair of neighbors
% %%
% % %     find the neighborhood at the image domain
% %     a = reshape(dist.flair_dist(mask),numel(dist.flair_dist(mask)), 1);
% %     IC = ones(size(mask))>0;
% %     inds = zeros(size(IC));
% %     iroi = find(IC);
% % %     iroi = a;
% %     inds(iroi) = dist.flair_dist;
% % %     inds(iroi) = 1:numel(iroi);
% %     
% %     neighbs26 = [10:18]; %14 is middle element
% % %     neighbs26 = [51:75];
% % 
% %     numNeighbs = numel(neighbs26);
% %     Neighbs = zeros(nnz(voxel_selection_mask), numNeighbs);
% %     
% %     for ii = 1:numel(neighbs26)
% %         ind = neighbs26(ii);
% % %         hh = zeros(repmat(5, [1,3]));
% %         hh = zeros(repmat(3, [1,3]));
% %         hh(ind) = 1;
% %         N_x = imfilter(inds, hh);
% %         Neighbs(:,ii) = N_x(voxel_selection_mask);
% %     end
% %     
% %     clear N_x;
% %     % replicative iroi boundary voxel processing
% %     tmp = repmat(a, [1, numel(neighbs26)]);
% %     Neighbs(Neighbs==0) = tmp(Neighbs==0);
% %     clear tmp;
     %%     load data into arrays for lrm
%     data = [ints3d.flair_s1ints(voxel_selection_mask), dist.flair_dist(voxel_selection_mask), ...
%         ints3d.t2_s1ints(voxel_selection_mask), dist.t2_dist(voxel_selection_mask), ...
%         ints3d.t1_s1ints(voxel_selection_mask),...
%         dist.t1_dist(voxel_selection_mask), Neighbs]; %6 central intensities plus texture properties
    data_norm = [dist.flair_dist(vsmask), data_norm];
    data_train  = [dist.flair_dist(vsmask_balanced), data_train];
%     data = Neighbs;
    size(data_norm)
%     gt = gt(voxel_selection_mask);
%     gt(gt == 2) = 0; % disregard improving lesions
    %% add to all data
    alltrainData(J) = {data_train};
    alltrainGroundTruth(J) = {labels_train};
    alltestData(J) = {data_norm};
    alltestGroundTruth(J) = {labels};
    %% save values for visualization/save_nii use later
    volumes(end+1) = {size(im)};
    masks(end+1) = {vsmask};
end 
    %% K-fold split
% folds = 3;
% indices = crossvalind('Kfold', patients, folds);
testAUCs = [];
trainingpredAll = []; % voxelwise predictions
testGroundTruthAll = [];
% folds = 5;
% indices = crossvalind('Kfold', pati ents, folds);
% indices = reshape(repmat([1,2,3], [5,1]), [15 1]);
 folds = 3;
    indices = reshape(repmat([1,2,3], [5,1]), [15 1]);
    indices = randsample(indices, length(indices));
for K = 1:folds
    testPatients = patients(indices == K);
    test=ismember(patients,testPatients);
    trainingPatients = patients(~test);
    %% init vars to add patient-by-patient
    trainingData = [];
    trainingGroundTruth = [];
    testData = [];
    testGroundTruth = [];
%     testData = {};
%     testGroundTruth = {};
    %% load data for training and test
    for I = 1:length(patients)
         if ismember(patients(I), trainingPatients)
            trainingData = [trainingData; alltrainData{I}];
            trainingGroundTruth = [trainingGroundTruth; alltrainGroundTruth{I}];
         else
             testData = [testData; alltestData{I}];
            testGroundTruth = [testGroundTruth; alltestGroundTruth{I}];
         end
    end
    %% svm classifier
%     svmmodel = fitcsvm(trainingData, trainingGroundTruth);
%     svm(end+1) = {svmmodel};
    %% elasticnet
%         trainingData(:,3:4) = [];
        [B, FitInfo] = lassoglm(double(trainingData),double(trainingGroundTruth),'binomial','CV',5,'Alpha',0.5);%
        lassoPlot(B,FitInfo,'PlotType','CV');
        
        index = FitInfo.IndexMinDeviance;%IndexMinMSE;
        wfeat = B(:,index);
        bar(wfeat)
        nnz(wfeat)
%         featureList(wfeat~=0)
    %% lasso training
    % tic
    % [B,FitInfo] = lasso(double(trainingData),double(trainingGroundTruth),'Alpha',0.5, 'CV',10);
    % lassoPlot(B,FitInfo,'PlotType','CV');
    % b = B(:,find(FitInfo.Lambda == FitInfo.LambdaMinMSE));
    % toc
%     %% run on test patients
%     legendstring = {};
%     %     testPatients = [1];
%     %     testPatients = [1,2,4,5,6,7,8,10,12,14,15,16,17,18,19];
%% elasticnet sanity check on training
%             XTest = trainingData;
%             indx = FitInfo.IndexMinDeviance;
%             B0 = B(:,indx);
%             nonzeros = sum(B0 ~= 0)
%             cnst = FitInfo.Intercept(indx);
%             weights = [cnst;B0];
%             tic()
%             trainingpred = glmval(weights,XTest,'logit'); %weights stored in B1
%             toc()
%             [X,Y,T,AUCROC,OPTROCPT] = perfcurve(double(trainingGroundTruth), trainingpred, 1,'xCrit','FPR','TVals',min(trainingpred):range(trainingpred)/100:max(trainingpred));%testGroundTruth{J}
%             AUCROC
            XTest = trainingData;
    XTest(isnan(XTest) |isinf(XTest))=0;
    indx = FitInfo.IndexMinDeviance;
    B0 = B(:,indx);
    nonzeros = sum(B0 ~= 0)
    cnst = FitInfo.Intercept(indx);
    B1 = [cnst;B0];
    tic()
    trainingpred = glmval(B1,XTest,'logit'); %weights stored in B1
    toc()
    [X,Y,T,AUCROC,OPTROCPT] = perfcurve(double(trainingGroundTruth), trainingpred, 1,'xCrit','FPR',...
        'TVals',min(trainingpred):range(trainingpred)/100:max(trainingpred));%testGroundTruth{J}
    fprintf('\nTraining AUC: %.2f\n',AUCROC)
    %% Testing voxelwise:
    XTest = testData;
    XTest(isnan(XTest) |isinf(XTest))=0;
    tic()
    testpred = glmval(B1,XTest,'logit'); %weights stored in B1
    toc()
        texturecoef(end+1) = {B1};

    [X,Y,T,AUCROC,OPTROCPT] = perfcurve(double(testGroundTruth), testpred, 1,'xCrit','FPR',...
        'TVals', min(testpred):range(testpred)/100:max(testpred) );%testGroundTruth{J}
    fprintf('\nTest (voxel-level) AUC: %.2f\n',AUCROC)
    trainingpredAll = [trainingpredAll;testpred]; 
    testGroundTruthAll = [testGroundTruthAll;testGroundTruth]; 
    %% TESTING patient-wise:
    for J = testPatients%
        XTest = alltestData{find(patients==J)};
        XTest(isnan(XTest) |isinf(XTest))=0;
        tic()
        testpred = glmval(B1,XTest,'logit'); %weights stored in B1
        toc()
        %     testpred = (data * b);
        % %         testpred = (data * b(2)) + b(1);
        %     testpred = 1./(1+exp(-testpred));
        %                     fprintf('end multiplication');
        
        %% create the volume
%         pred_vol = zeros(size(allMaskTest{J}));
%         pred_vol(allMaskTest{J}) = testpred;
        num = num2str(J);
         dummy_nii = load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_study1_T2Wreg.nii.gz'));
        pred_vol = zeros(size(dummy_nii.img));
%                 pred_vol(voxel_selection_mask>0) = testpred;
        mask = masks{find(patients == J)};
        pred_vol(mask>0) = testpred;

        %% save predictions as .nii
%         a=3
        temp = make_nii(pred_vol);    
        temp.hdr = dummy_nii.hdr;
        %temp.hdr.dime.bitpix = 256; % signed char
%         save_nii(temp, strcat('../MSpatientdata/patient', num, '/patient', num, '_9chtexturePred.nii.gz'));
            save_nii(temp, strcat('../genMSdata2/patient', num, '/patient', num, '_texturedflair_', num2str(iteration),'.nii.gz'));

        %%
        [X,Y,T,AUCROC,OPTROCPT] = perfcurve(double(alltestGroundTruth{find(patients == J)}), testpred, 1,'xCrit','FPR',...
            'TVals',min(testpred):range(testpred)/100:max(testpred));%testGroundTruth{J}
        fprintf('\nTest AUC: %.2f\n',AUCROC)
        testAUCs(J) = AUCROC;
    end
%% svm sanity check on training
%             trainingpred = predict(svmmodel, trainingData);
%             [X,Y,T,AUCROC,OPTROCPT] = perfcurve(double(trainingGroundTruth), trainingpred, 1,'xCrit','FPR','TVals',min(trainingpred):range(trainingpred)/100:max(trainingpred));%testGroundTruth{J}
%             AUCROC
%% run on test patients
% 
%     for J = 1:length(testPatients)
%         num = testPatients(J);
%         patient_number = num;
% %         s.patient_number = num;
%         num = num2str(num);
%         folder = strcat('../MSpatientdata/patient', num);
% % %         Get a list of all files in the folder with the desired file name pattern
% %         filePattern = fullfile(folder, '*.nii*'); % Change to whatever pattern you need.
% %         theFiles = dir(filePattern);
% %         for k = 1 : length(theFiles)
% %             baseFile = theFiles(k).name;
% %             fullFileName = fullfile(folder, baseFile);
% %             fprintf(1, 'Now reading %s\n', fullFileName);
% % %                 variable allocation
% %             if strfind(baseFile, '1_T1')
% %                 
% %                 s.t1_s1file = fullFileName;
% %             elseif strfind(baseFile, '1_T2')
% %                 s.t2_s1file = fullFileName;
% %             elseif strfind(baseFile, '1_FLAIR')
% %                 s.flair_s1file = fullFileName;
% %             elseif strfind(baseFile, '2_T1')
% %                 s.t1_s2file = fullFileName;
% %             elseif strfind(baseFile, '2_T2')
% %                 s.t2vs_s2file = fullFileName;
% %             elseif strfind(baseFile, '2_FLAIR')
% %                 s.flair_s2file = fullFileName;
% %             elseif strfind(baseFile, 'gt3')
% %                 s.gtfile = fullFileName;
% %             elseif strfind(baseFile, 'mask')
% %                 s.maskfile = fullFileName;
% %             end
% %         end
% %         % get test data on the whole brainmask
% % 
% %         mask = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_brainmask.nii.gz')),'img')>0;
% % 
% %         gt = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient', num,'_gt3.nii')),'img');
% %         gt(gt==2) = 0;%    
% %          fields = fieldnames(s);
% %     fields = setdiff(fields,{'maskfile','gtfile','patient_number'});%%skip patient number, also skip the brain mask and the ground truth mask
% %     ints3d = struct;
% % 
% %     if strcmp(normmethod,'zscore')
% %         for N = 1:numel(fields)
% %             field = fields{N};
% %             value = getfield(s, field);
% %             nii = load_nii(value);
% %             im = nii.img;
% %             r = strrep(field,'file','');
% %             temp = im(mask);
% %             centered = (temp - mean(temp)) ;
% %             %         sum = 0;
% %             %         for i = 1:numel(centered)
% %             %             sum = sum + (centered(i, 1)^2);
% %             %         end
% %             stddev = std(double(temp));
% %             %     mean(temp)
% %             %     ints.(strcat(r,'ints')) = centered; %0 centered, now all stored in ints - should this be before or after difference calculation?
% %             ints3d.(strcat(r,'ints')) = zeros(size(mask));
% %             ints3d.(strcat(r,'ints'))(mask) = centered/stddev;
% %         end
% %        %% normalize the intensities wrt WM
% %     elseif strcmp(normmethod,'wm')
% %         [ints3d.t1_s1ints, ints3d.t2_s1ints, ints3d.flair_s1ints] = normalizeWM(getfield(load_nii(s.t1_s1file),'img'),...
% %             getfield(load_nii(s.t2_s1file),'img'),...
% %             getfield(load_nii(s.flair_s1file),'img'), mask);
% %     end
%    
%     %% calc distances
% %     dist.t1_dist = ints3d.t1_s2ints - ints3d.t1_s1ints;
% %     dist.flair_dist = ints3d.flair_s2ints - ints3d.flair_s1ints;
% %     dist.t2_dist = ints3d.t2_s2ints - ints3d.t2_s1ints;
% % %      %% 
% %      % find the neighborhood at the image domain
% %         a = reshape(dist.flair_dist(mask),numel(dist.flair_dist(mask)), 1);
% %     IC = ones(size(mask))>0;
% %     inds = zeros(size(IC));
% %     iroi = find(IC);
% %     %     iroi = a;
% %     inds(iroi) = dist.flair_dist;
% %     %     inds(iroi) = 1:numel(iroi);
% % 
% % %     neighbs26 = [10,11,12,13,15,16,17,18];
% %         neighbs26 = [10:18];
% % %     neighbs26 = [51:75];
% %     numNeighbs = numel(neighbs26);
% %     Neighbs = zeros(nnz(mask), numNeighbs);
% % 
% %     for ii = 1:numel(neighbs26)
% %         ind = neighbs26(ii);
% %         hh = zeros(repmat(3, [1,3]));
% %         hh(ind) = 1;
% %         N_x = imfilter(inds, hh);
% %         Neighbs(:,ii) = N_x(mask);
% %     end
% % 
% %     clear N_x;  
%     %% compute texture features
% %                 gt = getfield(load_nii([folder,'/patient', num,'_gt3.nii']),'img');
% 
% %      voxel_selection_mask = zeros(size(gt));
% %     voxel_selection_mask(gt==1) = 1; % select all positive class observations
% %     voxel_selection_mask(randsample(find(mask>0 & gt~=1 & ints3d.t2_s1ints>0),nnz(gt==1))) = 1; % randomly select the same number of other observations
% %     voxel_selection_mask = logical(voxel_selection_mask);   
% %     textureProps =[];
% %     patchSize = 21;
% %     empty = ones(patchSize);
% % %     empty((patchSize * patchSize + 1)/2) = 1;
% %     for i=1:size(dist.flair_dist,3) %by slice
% %             patches = nlfilter(dist.flair_dist(:,:,i), [patchSize patchSize], @(block) {block});
% % %             maskslice = voxel_selection_mask(:,:,i);
% %             maskslice = mask(:,:,i);
% %             patches(maskslice == 0) = []; 
% %             for k=1:numel(patches)
% %                 temp = patches{k};
% %                 [feat,featureGroupID, featureGroupList, featureList] = cellfun(@computeFeats, patches, 'UniformOutput', false);
% %                 textureProps = [textureProps, [feat{:}]];
% %                 glcm = graycomatrix(temp);
% %                 grayprops = graycoprops(glcm);
% % %                 temp = reshape(patches(:,k), patchSize, patchSize);
% %                 [feat,featureGroupID, featureGroupList, featureList] = computeFeats(temp, logical(empty));
%       
% %             end
% %     end   
% %   
% %                 [feat,featureGroupID, featureGroupList, featureList] = computeFeats(temp, empty==1);
% %                 textureProps = [textureProps, feat];
% %         textureProps = [textureProps; grayprops.Contrast,grayprops.Correlation,grayprops.Energy,grayprops.Homogeneity]; %add graycoprops value(s) to data matrix
% 
%             
% %     end
%     %% fix texture properties: elim NaN and normalize
% %     textureProps(isnan(textureProps))=0;
% %     textureProps = textureProps';
% %     %% normalize scaling between 0 and 1
% % %     denominator = max(textureProps) - min(textureProps);
% % %     sub = min(textureProps);
% % %% normalize using z-scores
% %     denominator = std(textureProps);
% %     sub = mean(textureProps);
% %     %% appply normalization
% %     denominator = repmat(denominator, [size(textureProps,1),1]);
% %     sub = repmat(sub, [size(textureProps,1),1]);
% %     result = (textureProps - sub)./denominator;
% % 
% %             %%             load data into arrays for lrm
% % %             data = [ints3d.flair_s1ints(mask), dist.flair_dist(mask), ints3d.t2_s1ints(mask), dist.t2_dist(mask), ints3d.t1_s1ints(mask), dist.t1_dist(mask), Neighbs];
% % %     data = Neighbs;
% %     data = [dist.flair_dist(mask), result];
% %         data(:,3:4) = [];
% 
%         %% get  predictions
%         
% %         a=1
%     %     patient_number = 1;
%         %patient_number = testPatients(J); %get patient's number
%         data = testData{J};
%                 fprintf('start multiplication');
%                 %% get elasticnet predictions
% %         indstest = find(Indices==ifold);
% %             XTest =  FEATS(indstest,:);
%             XTest = data;
%             indx = FitInfo.IndexMinDeviance;
%             B0 = B(:,indx);
%             nonzeros = sum(B0 ~= 0)
%             cnst = FitInfo.Intercept(indx);
%             B1 = [cnst;B0];
%             tic()
%             testpred = glmval(B1,XTest,'logit'); %weights stored in B1
%             toc()
% %% get lasso/glmfit predictions
% 
%     %     testpred = (data * b);
%     % %         testpred = (data * b(2)) + b(1);
%     %     testpred = 1./(1+exp(-testpred));
%     %                     fprintf('end multiplication');
%     %% get svm predictions
% %     testpred = predict(svmmodel, data);
%    
%         %% create the volume
    %     testpred3d = ind2sub(testVolume(J),testpred);
%     a=2
%         dummy_nii = load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_study1_T2Wreg.nii.gz'));
%         pred_vol = zeros(size(dummy_nii.img));
% %                 pred_vol(voxel_selection_mask>0) = testpred;
%         mask = masks{find(patients == patient_number)};
%         pred_vol(mask>0) = testpred;
% 
%         %% save predictions as .nii
% %         a=3
%         temp = make_nii(pred_vol);    
%         temp.hdr = dummy_nii.hdr;
%         %temp.hdr.dime.bitpix = 256; % signed char
%         save_nii(temp, strcat('../MSpatientdata/patient', num, '/patient', num, '_dftexturePred.nii.gz'));
%         %% calculate AUC 
% % a=4
%     %         % plot the ROC curve:
%     %         subplot(1,2,1)
%     %         [X,Y,T,AUCROC,OPTROCPT] = perfcurve(gt(mask>0), testpred, 1,'xCrit','FPR');%testGroundTruth{J}
%     %         AUCROC
%     %         testAUCROC(s.patient_number) = AUCROC;
%     %         hold on;
%     %         legendstring{end+1} = strcat('Patient ', num);
%     %         plot(X,Y);
%     %         legend(legendstring);
%     %     %     title = strcat
%     %         title('Test Data ROC');
%     %         axis square
% 
%     %         % plot the ROC curve:
%     %         subplot(1,2,2)
%     %         [X,Y,T,AUCPR,OPTROCPT] = perfcurve(gt(mask>0), testpred, 1,'xCrit','PPV');%testGroundTruth{J
%     %         AUCPR
%     %         testAUCPR(s.patient_number) = AUCPR;
%     %         hold on;
%     %         legendstring{end+1} = strcat('Patient ', num);
%     %         plot(X,Y);
%     %         legend(legendstring);
%     %     %     title = strcat
%     %         title('Test Data PR');
%     %         axis square
%     %         %% dsi
%     %         dice = @(auto_seg,manual_seg) 2*nnz(auto_seg & manual_seg)/(nnz(auto_seg) + nnz(manual_seg));
%     %         dsitest = testpred;
%     %         dsitest(dsitest > 0.5) = 1;
%     %         dsitest(dsitest < 0.5) = 0;
%     %         dsi = dice(dsitest, gt(mask>0));
%     %         
%     %         dsi
%     %         testDSI(s.patient_number) = dsi;
%     end
% %     B1
%     texturecoef(end+1) = {B1};

end
%% Overall voxelwise performance:
[X,Y,T,AUCROC,OPTROCPT] = perfcurve(double(testGroundTruthAll), trainingpredAll, 1,'xCrit','FPR',...
    'TVals', min(trainingpredAll):range(trainingpredAll)/100:max(trainingpredAll) );%testGroundTruth{J}
fprintf('\nTest (voxel-level) AUC: %.2f\n',AUCROC)
%%
fprintf('\nTest (patient-level) AUC mean: %.2f(%.2f)\n',mean(testAUCs(patients)),std(testAUCs(patients)))
end
%% post processing
% imgs = {'dect', 'manual', 'dFlair'};
imgs = {'9chtexturePred'}
% patients = [2]
% T = struct;
for i = 1:length(imgs)
    name = imgs{i}
    auc = [];
    opt = [];
    for J = 1:length(patients)
        %% load files
        patient_number = patients(J); %change the particular patient's number
        num = num2str(patient_number);
        mask = masks{find(patients == patient_number)};

%         mask = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_brainmask.nii.gz')),'img')>0;
        gt3 = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_gt3.nii')),'img')==1;
        img = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_',name,'.nii.gz')),'img');
        %% voxel selection mask
%         voxel_selection_mask = zeros(size(gt3));
%         voxel_selection_mask(img~=0) = 1; % select all positive class observations
%         voxel_selection_mask(find(mask>0 & gt3~=1 & img>0)) = 1; 
%         voxel_selection_mask(find(mask>0 & gt3~=1 & img>0),nnz(gt3==1))) = 1; 
%         voxel_selection_mask = logical(voxel_selection_mask);
        %% 
%         new = img(voxel_selection_mask);
        new = img(mask);
%         gt = gt3(voxel_selection_mask);
        gt = gt3(mask);
        new(isnan(new) |isinf(new))=0;

        [X,Y,T,AUCROC,OPTROCPT] = perfcurve(gt, new, true,'xCrit','FPR','TVals',min(new):range(new)/100:max(new));%testGroundTruth{J}
        AUCROC;
        auc = [auc; AUCROC];
        thr = T(find(X>=OPTROCPT(1),1));
        opt = [opt; thr];
        plot(X,Y);
        axis square
    end
    auc
    opt
    if strfind(name, 'manual')
        T.manual_dect_auc = auc';
        T.manual_dect_opt = opt';
    elseif strfind(name, 'dect')
        T.trained_dect_auc = auc';
        T.trained_dect_opt = opt';
    elseif strfind(name, 'Flair')
        T.flair_auc = auc';
        T.flair_opt = opt';
    end
    
end
%% if filename wrong
%     for J = 1:length(patients)
%         patient_number = patients(J); %change the particular patient's number
%         num = num2str(patient_number);
%         img = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_','texturePred.nii.gz')),'img');
% 
%           dummy_nii = load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_texturePred.nii.gz'));
% 
%         %% save predictions as .nii
% %         a=3
%         temp = make_nii(img);    
%         temp.hdr = dummy_nii.hdr;
%         %temp.hdr.dime.bitpix = 256; % signed char
%         save_nii(temp, strcat('../MSpatientdata/patient', num, '/patient', num, '_dFlairNeighborPred.nii.gz'));
% 
%     end

