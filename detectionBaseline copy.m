for iteration = 1:3
    % sweeney-based LRM change detection baseline
    %% init vars
    bverbose = true;
    normmethod = 'zscore'; % normalization method: 'zscore','wm'
    testAUC = [];
    allData = {};
    volumes = {};
    masks = {};
    allGroundTruth = {};
    patients = [1,2,4,5,6,7,8,10,12,14,15,16,17,18,19]; %no AUC = 0 no disease progression from combine.m
    % patients = [1,2,4,5,6,7,8,10,12,14,15,16,17,18,19]; %no AUC = 0 no disease progression from combine.m
    % trainingPatients = [1,2,4,5,6,7,8,10,12,14];
    % testPatients = [15,16,17,18,19];
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
        elseif strcmp(normmethod,'wm')
            %% normalize the intensities wrt WM
            [ints3d.t1_s1ints, ints3d.t2_s1ints, ints3d.flair_s1ints] = normalizeWM(getfield(load_nii(s.t1_s1file),'img'),...
                getfield(load_nii(s.t2_s1file),'img'),...
                getfield(load_nii(s.flair_s1file),'img'), mask);
        end
        %% calc distances

        dist.t1_dist = ints3d.t1_s2ints - ints3d.t1_s1ints;
        dist.flair_dist = ints3d.flair_s2ints - ints3d.flair_s1ints;
        dist.t2_dist = ints3d.t2_s2ints - ints3d.t2_s1ints;
        %%
    %     voxel_selection_mask = zeros(size(ints3d.t2_s1ints));
    %     voxel_selection_mask(ints3d.t2_s1ints>0) =1;
    %     voxel_selection_mask = logical(voxel_selection_mask);
    %% new dflair vsm
        temp = zeros(size(dist.flair_dist));
        vsmask = zeros(size(dist.flair_dist));
        for i=1:size(dist.flair_dist,3)
            temp(:,:,i) = medfilt2(dist.flair_dist(:,:,i));
        end
        vsmask(temp>(std(temp(:)))) = 1;

        %%
        voxel_selection_mask = zeros(size(gt));
        voxel_selection_mask(gt==1) = 1; % select all positive class observations
        voxel_selection_mask(randsample(find(mask>0 & gt~=1 &vsmask==1),nnz(gt==1))) = 1; % randomly select the same number of other observations
        voxel_selection_mask = logical(voxel_selection_mask);
        %%     load data into arrays for lrm
    %     data = [dist.flair_dist(voxel_selection_mask)];
        data = [ints3d.flair_s1ints(voxel_selection_mask), dist.flair_dist(voxel_selection_mask), ints3d.t2_s1ints(voxel_selection_mask), dist.t2_dist(voxel_selection_mask), ints3d.t1_s1ints(voxel_selection_mask), dist.t1_dist(voxel_selection_mask)];
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
        %% K-fold split
    folds = 3;
    indices = reshape(repmat([1,2,3], [5,1]), [15 1]);
    indices = randsample(indices, length(indices));
%     indices = crossvalind('Kfold', patients, folds);
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
             if ismember(patients(I), trainingPatients)
                trainingData = [trainingData; allData{I}];
                trainingGroundTruth = [trainingGroundTruth; allGroundTruth{I}];
            else
        %         testData = [testData; data];
        %         testGroundTruth = [testGroundTruth; gt];
                testData(end+1) = {allData{I}};
                testGroundTruth(end+1) = {allGroundTruth{I}};
             end
        end
        %% lasso training
        % tic
        % [B,FitInfo] = lasso(double(trainingData),double(trainingGroundTruth),'CV',10);
        % lassoPlot(B,FitInfo,'PlotType','CV');
        % b = B(:,find(FitInfo.Lambda == FitInfo.LambdaMinMSE, size(trainingData,2)));
        % toc
        %% elasticnet
    %     [B, FitInfo] = lassoglm(double(trainingData),double(trainingGroundTruth),'binomial','CV',5,'Alpha',0.5);%
    %             % lassoPlot(B,FitInfo,'PlotType','CV');
    %             %
    %             index = FitInfo.IndexMinDeviance;%IndexMinMSE;
    %             wfeat = B(:,index);
    %             bar(wfeat)
    %             nnz(wfeat)
        %% glmfit training
            tic
            [b, dev, stats] = glmfit(double(trainingData),double(trainingGroundTruth), 'normal');
            fprintf('Glmfit ');
            toc

        %% create predictions for each training patient
        %     training
        %     trainingpred = (trainingData * b(2:7)) + b(1);
        % new = zeros(size(mask));
        %     trainingpred = 1./(1+exp(-trainingpred));
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
        %     testPatients = [1];
    %     testPatients = [1,2,4,5,6,7,8,10,12,14,15,16,17,18,19];
        for J = 1:length(testPatients)
            num = testPatients(J);
            s.patient_number = num;
            num = num2str(num);
            folder = strcat('../MSpatientdata/patient', num);
            % Get a list of all files in the folder with the desired file name pattern
            filePattern = fullfile(folder, '*.nii*'); % Change to whatever pattern you need.
            theFiles = dir(filePattern);
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

            elseif strcmp(normmethod,'wm')
                %% normalize the intensities wrt WM
                [ints3d.t1_s1ints, ints3d.t2_s1ints, ints3d.flair_s1ints] = normalizeWM(getfield(load_nii(s.t1_s1file),'img'),...
                    getfield(load_nii(s.t2_s1file),'img'),...
                    getfield(load_nii(s.flair_s1file),'img'), mask);
            end
            %% calc distances
            dist.t1_dist = ints3d.t1_s2ints - ints3d.t1_s1ints;
            dist.flair_dist = ints3d.flair_s2ints - ints3d.flair_s1ints;
            dist.t2_dist = ints3d.t2_s2ints - ints3d.t2_s1ints;
            %%

    %             voxel_selection_mask = zeros(size(ints3d.t2_s1ints));
    %             voxel_selection_mask(ints3d.t2_s1ints>0) =1;
        %         voxel_selection_mask = logical(voxel_selection_mask);
        %         mask = voxel_selection_mask;
            %%     load data into arrays for lrm
    %           data = dist.flair_dist(mask);
                data = [ints3d.flair_s1ints(mask), dist.flair_dist(mask), ints3d.t2_s1ints(mask), dist.t2_dist(mask), ints3d.t1_s1ints(mask), dist.t1_dist(mask)];
        %         b = [-9.1008; 0.7388; -0.054; -0.2531; 0.6503; 0.5098; -0.8282];
            %% get predictions
        %     patient_number = 1;
            %patient_number = testPatients(J); %get patient's number
    %         testpred = data*b;
                testpred = (data * b(2:end)) + b(1);
    %             testpred = (data * b(2)) + b(1);

            testpred = 1./(1+exp(-testpred));
            %% get predictions for elasticnet
    %           XTest = data;
    %                 indx = FitInfo.IndexMinDeviance;
    %                 B0 = B(:,indx);
    %                 nonzeros = sum(B0 ~= 0)
    %                 cnst = FitInfo.Intercept(indx);
    %                 B1 = [cnst;B0];
    %                 tic()
    %                 testpred = glmval(B1,XTest,'logit');
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
            save_nii(temp, strcat('../genMSdata2/patient', num, '/patient', num, '_baselinet_', num2str(iteration),'.nii.gz'));
            %% calculate AUC 
            gt = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient', num,'_gt3.nii')),'img');
            gt(gt==2) = 0;
        %         % plot the ROC curve:
        %         subplot(1,2,1)
        %         [X,Y,T,AUCROC,OPTROCPT] = perfcurve(gt(mask>0), testpred, 1,'xCrit','FPR');%testGroundTruth{J}
        %         AUCROC
        %         testAUCROC(s.patient_number) = AUCROC;
        %         hold on;
        %         legendstring{end+1} = strcat('Patient ', num);
        %         plot(X,Y);
        %         legend(legendstring);
        %     %     title = strcat
        %         title('Test Data ROC');
        %         axis square

        %         % plot the ROC curve:
        %         subplot(1,2,2)
        %         [X,Y,T,AUCPR,OPTROCPT] = perfcurve(gt(mask>0), testpred, 1,'xCrit','PPV');%testGroundTruth{J
        %         AUCPR
        %         testAUCPR(s.patient_number) = AUCPR;
        %         hold on;
        %         legendstring{end+1} = strcat('Patient ', num);
        %         plot(X,Y);
        %         legend(legendstring);
        %     %     title = strcat
        %         title('Test Data PR');
        %         axis square
            %% dsi
        %         dice = @(auto_seg,manual_seg) 2*nnz(auto_seg & manual_seg)/(nnz(auto_seg) + nnz(manual_seg));
        %         dsitest = testpred;
        %         dsitest(dsitest > 0.5) = 1;
        %         dsitest(dsitest < 0.5) = 0;
        %         dsi = dice(dsitest, gt(mask>0));
        %         
        %         dsi
        %         testDSI(s.patient_number) = dsi;
        end
    end
end
%% post processing
imgs = {'vsmask'}
% imgs = {'dect', 'manual', 'dFlair'};
% T = struct;
for i = 1:length(imgs)
    name = imgs{i}
    auc = [];
    opt = [];
    for J = 1:length(patients)
        %% load files
        patient_number = patients(J); %change the particular patient's number
        num = num2str(patient_number);
        mask = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_brainmask.nii.gz')),'img')>0;
        gt3 = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_gt3.nii')),'img')==1;
        img = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_',name,'Baselinetest.nii.gz')),'img');
        %% voxel selection mask
        voxel_selection_mask = zeros(size(gt3));
        voxel_selection_mask(gt3==1) = 1; % select all positive class observations
        voxel_selection_mask(randsample(find(mask>0 & gt3~=1 & img>0),nnz(gt3==1))) = 1; 
        voxel_selection_mask = logical(voxel_selection_mask);
        %% 
        new = img(voxel_selection_mask);
        gt = gt3(voxel_selection_mask);
        [X,Y,T,AUCROC,OPTROCPT] = perfcurve(gt, new, true,'xCrit','FPR','TVals',min(new):range(new)/100:max(new));%testGroundTruth{J}
        AUCROC
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

%% save table
% writetable(T,'baselineAUC.csv','Delimiter',',','QuoteStrings',true)
% type 'baselineAUC.csv'

        %



% % %%
% % [prec, tpr, fpr, thresh] = prec_rec( gt(mask>0)*0.5, gt(mask>0)>0 );%testpred
% % prec_rec( double(gt(mask>0))*0.5, gt(mask>0)>0 );%testpred
% % legendstring{end+1} = strcat('Patient ', num);
% % legend(legendstring);
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

