% voxel-level prediction baseline
%% init vars
bverbose = true;
normmethod = 'zscore'; % normalization method: 'zscore','wm'
testAUC = [];
allData = {};
volumes = {};
masks = {};
allGroundTruth = {};
patients = [1,2,4,5,6,7,8,10,12,14,15,16,17,18,19]; %no AUC = 0 no disease progression from combine.m
%% manual holdout
% testPatients = datasample(patients,3,'Replace',false);
% test=ismember(patients,testPatients);
% trainingPatients = patients(~test);

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
        elseif strfind(baseFile, 'gt3')
            s.gtfile = fullFileName;
        elseif strfind(baseFile, 'mask')
            s.maskfile = fullFileName;
        end
        %     nii = load_nii(fullFileName);
        %     if bview
        %         view_nii(nii);
        %         drawnow; % Force display to update immediately.
        %     end
    end
%     s
    %% load and normalize data

    nii = load_nii(s.maskfile);
    mask = nii.img;
    mask = logical(mask);
    nii = load_nii(s.gtfile);
    %gt = nii.img; % can do this: getfield(load_nii(s.gtfile),'img')
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
        % dist.t1_dist = ints.t1_s2ints - ints.t1_s1ints;
        % dist.t2_dist = ints.t2_s2ints - ints.t2_s1ints;
    elseif strcmp(normmethod,'wm')
        %% normalize the intensities wrt WM
        [ints3d.t1_s1ints, ints3d.t2_s1ints, ints3d.flair_s1ints] = normalizeWM(getfield(load_nii(s.t1_s1file),'img'),...
            getfield(load_nii(s.t2_s1file),'img'),...
            getfield(load_nii(s.flair_s1file),'img'), mask);
    end
    %%
%     voxel_selection_mask = zeros(size(ints3d.t2_s1ints));
%     voxel_selection_mask(ints3d.t2_s1ints>0) =1;
%     voxel_selection_mask = logical(voxel_selection_mask);
    %%
    voxel_selection_mask = zeros(size(gt));
    voxel_selection_mask(gt==1) = 1; % select all positive class observations
    voxel_selection_mask(randsample(find(mask>0 & gt~=1 & ints3d.t2_s1ints>0),nnz(gt==1))) = 1; % randomly select the same number of other observations
    voxel_selection_mask = logical(voxel_selection_mask);
    %%     load data into arrays for lrm
    data = [ints3d.flair_s1ints(voxel_selection_mask), ints3d.t2_s1ints(voxel_selection_mask), ints3d.t1_s1ints(voxel_selection_mask)];
    gt = gt(voxel_selection_mask);
    gt(gt == 2) = 0; % disregard improving lesions
    allData(J) = {data};
    allGroundTruth(J) = {gt};
    %% add to array
%     allData = [allData; data];
%     allGroundTruth = [allGroundTruth; gt];
    %% save values for visualization/.nii save later
    volumes(end+1) = {size(im)};
    masks(end+1) = {mask};
    %% patient level holdout set up input data for training by adding to input matrix (truths = gt(mask))

%     if ismember(patient_number, trainingPatients)
%         trainingData = [trainingData; data];
%         trainingGroundTruth = [trainingGroundTruth; gt];
%     else
% %         testData = [testData; data];
% %         testGroundTruth = [testGroundTruth; gt];
%         testVolume(end+1) = {size(im)};
%         testData(end+1) = {data};
%         testGroundTruth(end+1) = {gt};
%         testMask(end+1) = {mask};
% 
%     end
end
%% %% K-fold split
folds = 3;
indices = crossvalind('Kfold', patients, folds);
%% run model using K-fold
figure;

for K = 1:folds
    testPatients = patients(indices == K);
    test=ismember(patients,testPatients);
    trainingPatients = patients(~test);
    %% init vars to add patient-by-patient
    trainingData = [];
    trainingGroundTruth = [];
    testData = {};
    testGroundTruth = {};
    %% load data for training and test
    for I = 1:length(patients)
         if ismember(I, trainingPatients)
            trainingData = [trainingData; allData{I}];
            trainingGroundTruth = [trainingGroundTruth; allGroundTruth{I}];
        else
    %         testData = [testData; data];
    %         testGroundTruth = [testGroundTruth; gt];
            testData(end+1) = {allData{I}};
            testGroundTruth(end+1) = {allGroundTruth{I}};
         end
    end
    %% training
    tic
    [b, dev, stats] = glmfit(double(trainingData),double(trainingGroundTruth), 'normal');
    fprintf('Glmfit ');
    toc
    % %% glm val
    % yfit = glmval(b,x,'probit','size',n);
    % plot(x, y./n,'o',x,yfit./n,'-','LineWidth',2)

    %% create predictions for each training patient
    trainingpred = (trainingData * b(2:4)) + b(1);
    % new = zeros(size(mask));
    trainingpred = 1./(1+exp(-trainingpred));
    % trainingpred(trainingpred<0.5) = 0;
    % trainingpred(trainingpred>0.5) = 1;
    % new(mask) = trainingpred';
    %% AUC and ROC for training results
    % tic
    % [X,Y,T,trainingAUC] = perfcurve(trainingGroundTruth, trainingpred, 1);
    % trainingAUC
    % figure;
    % plot(X,Y);
    % title('Training Data ROC');
    % fprintf('TrainingAUC ');
    % toc
    %% visualize
    % islice = 38;
    % % subplot(1,2,1);
    % imagesc(new(:,:,islice), [0, 1]);
    % axis image;
    %% run on test patients
    legendstring = {};
    for J = 1:length(testPatients)
        
        num = testPatients(J};
        num = num2str(num);
        s.patient_number = num;
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
            elseif strfind(baseFile, 'gt3')
                s.gtfile = fullFileName;
            elseif strfind(baseFile, 'mask')
                s.maskfile = fullFileName;
            end
        end
        %% get test data on the whole brainmask
        
        mask = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_brainmask.nii.gz')),'img')>0;
%         ints3d.t1_s1ints = load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_study1_T1Wreg.nii.gz'));
%         ints3d.t2_s1ints = load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_study1_T2Wreg.nii.gz'));
%         ints3d.flair_s1ints = load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_study1_FLAIRWreg.nii.gz'));
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
            % dist.t1_dist = ints.t1_s2ints - ints.t1_s1ints;
            % dist.t2_dist = ints.t2_s2ints - ints.t2_s1ints;
        elseif strcmp(normmethod,'wm')
            %% normalize the intensities wrt WM
            [ints3d.t1_s1ints, ints3d.t2_s1ints, ints3d.flair_s1ints] = normalizeWM(getfield(load_nii(s.t1_s1file),'img'),...
                getfield(load_nii(s.t2_s1file),'img'),...
                getfield(load_nii(s.flair_s1file),'img'), mask);
        end
        %%
%         voxel_selection_mask = zeros(size(ints3d.t2_s1ints));
%         voxel_selection_mask(ints3d.t2_s1ints>0) =1;
%         voxel_selection_mask = logical(voxel_selection_mask);
%         mask = voxel_selection_mask;
        %%     load data into arrays for lrm
        data = [ints3d.flair_s1ints(mask), ints3d.t2_s1ints(mask), ints3d.t1_s1ints(mask)];
        
        %% get predictions
    %     patient_number = 1;
        %patient_number = testPatients(J); %get patient's number
        testpred = (data * b(2:4)) + b(1);
        testpred = 1./(1+exp(-testpred));
        %% visualize
%         new = zeros(volumes{patient_number});
%         new(masks{patient_number}) = testpred;
%         size(new)
%     %     testpred3d = ind2sub(testVolume{J},testpred);
%     % set colors to define each one
%         islice = 27;
%     % % subplot(1,8,1);
%         imagesc(new(:,:,islice), [0, 1]);
%         axis image;
        %% dummy predictor
%         testpred = ints3d.t2_s1ints(mask);
%         cands = gt(mask)==1;%ints3d.t2_s1ints(mask)>0;
%         testpred(~cands) = 0;
%         testpred(cands) = 0.9;%rand([nnz(cands),1]);
%         inds = find(~cands);
%         testpred(inds(1)) = 0.9;
        %% create the volume
    %     testpred3d = ind2sub(testVolume(J),testpred);
        dummy_nii = load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_study1_T2Wreg.nii.gz'));
        pred_vol = zeros(size(dummy_nii.img));
        pred_vol(mask>0) = testpred;
        
        %% save predictions as .nii
        temp = make_nii(pred_vol);    
        temp.hdr = dummy_nii.hdr;
        %temp.hdr.dime.bitpix = 256; % signed char
        save_nii(temp, strcat('../MSpatientdata/patient', num, '/patient', num, '_predBaselinetest.nii.gz'));
        %% calculate AUC 
        gt = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient', num,'_gt3.nii')),'img');
        gt(gt==2) = 0;
        % plot the ROC curve:
        subplot(1,2,1)
         [X,Y,T,AUCROC,OPTROCPT] = perfcurve(gt(mask>0), testpred, 1,'xCrit','FPR');%testGroundTruth{J}
        AUCROC
        testAUCROC(s.patient_number) = AUCROC;
        hold on;
        legendstring{end+1} = strcat('Patient ', num);
        plot(X,Y);
        legend(legendstring);
    %     title = strcat
        title('Test Data ROC');
        axis square
        
        % plot the ROC curve:
        subplot(1,2,2)
        [X,Y,T,AUCPR,OPTROCPT] = perfcurve(gt(mask>0), testpred, 1,'xCrit','PPV');%testGroundTruth{J
        AUCPR
        testAUCPR(s.patient_number) = AUCPR;
        hold on;
        legendstring{end+1} = strcat('Patient ', num);
        plot(X,Y);
        legend(legendstring);
    %     title = strcat
        title('Test Data PR');
        axis square
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the delta fl
Xaxis = [];
Yaxis = [];
for J = 1:length(patients)
    num = num2str(patients(J));
    mask = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_brainmask.nii.gz')),'img')>0;
        gt3 = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_gt3.nii')),'img')==1;
        img = getfield(load_nii(strcat('../MSproject/zscore/patient',num,'/patient',num,'_delta_FLAIR.nii.gz')),'img');
        img2 = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient', num,'_dFlairBaselinetest.nii.gz')),'img');
        Xaxis = [Xaxis; img(mask)];
                Yaxis = [Yaxis; img2(mask)];
end
plot(Xaxis, Yaxis);
title('Delta_flair vs predictions on delta_flair across all patients');
%% mask for AUC
    new = img(mask);
    gt = gt3(mask);
% %%
% voxel_selection_mask = zeros(size(gt3));
% voxel_selection_mask(gt3==1) = 1; % select all positive class observations
% voxel_selection_mask(randsample(find(mask>0 & gt3~=1 & img>0),nnz(gt3==1))) = 1; 
% voxel_selection_mask = logical(voxel_selection_mask);
% %%
% new = img(voxel_selection_mask);
% gt = gt3(voxel_selection_mask);
%%
[X,Y,T,AUCROC,OPTROCPT] = perfcurve(gt, new, true,'xCrit','PPV','TVals',0:0.01:3);%testGroundTruth{J}
AUCROC
plot(X,Y);
axis square
%%
% homemade roc curve for validation only:
figure
a  = 'hi'
xfun = @(thr) 1-(nnz(~(gt > 0) & new<thr)/nnz(~(gt > 0))); %1-sp
xfun2 = @(thr) ( (nnz(gt & new>=thr)) / (nnz(new>=thr)) ); %ppv
yfun = @(thr) nnz(gt & new>=thr)/nnz(gt); %se
x = []; y=[];
T = [];
for thr = 0:0.01:3%min(new):range(new)/100:max(new) %unique(new(:)')%[2,3,4,7]%unique(TIRADS)'%
    x = [x;xfun2(thr)];
    y = [y;yfun(thr)];
    T = [T;thr];
end
% x = [x;1];y = [y;0];
%[~,ix] = sort(y);
clf, plot(x(:),y(:),'.-')
axis square
AUC = 0.5*sum( (x(2:end)-x(1:end-1)).*(y(2:end)+y(1:end-1)) )
% AUC = abs(AUC);

%%  
% [prec, tpr, fpr, thresh] = prec_rec( new, gt>0 );%testpred
prec_rec( new, gt>0 );%testpred
legendstring{end+1} = strcat('Patient ', num);
legend(legendstring);
%%
% testpred = (testData * b(2:4)) + b(1);
% new = zeros(size(mask));
% testpred = 1./(1+exp(-testpred));
% testpred(testpred<0.5) = 0;
% testpred(testpred>0.5) = 1;
% new(mask) = testpred';
%% validate performance
% % preds{J} = new;
% [X,Y,T,AUC] = perfcurve(testGroundTruth, testpred, 1);
% %new(mask) just shows whether it is part of the brain mask or not
% AUC    % auc(J) = AUC;
% figure;
% plot(X,Y)
% title('Test Data ROC');

