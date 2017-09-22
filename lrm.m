%%TODO: save and visualize .nii, pr-roc,
addpath(genpath('./deps'))
%% create the new (three-label) gts (NB: needed only once!):
for J = 1:20
    [gt, mask] = combine(J, false, true);
end
%% create the wm probability masks:

%%
tic
preds = cell(20, 1);
gts = cell(20, 1);
rocauc = zeros(20, 1);
prauc = zeros(20, 1);

bverbose = true;
normmethod = 'zscore'; % normalization method: 'zscore','wm'
patients = [1]
% patients = [1,2,4,5,6,7,8,10,12,14,15,16,17,18,19]; %no AUC = 0 no disease progression from combine.m

for J = 1:length(patients)
    %% loading data into variables for a patient to create structure
    patient_number =  patients(J); %change the particular patient's number
    num = num2str(patient_number);
    folder = strcat('../MSpatientdata/patient', num);
    % Get a list of all files in the folder with the desired file name pattern
    filePattern = fullfile(folder, '*.nii.gz'); % Change to whatever pattern you need.
    theFiles = dir(filePattern);
    % data = ones(2, 3);
    s.patient_number = num;
    
    for k = 1 : length(theFiles)
        baseFile = theFiles(k).name;
        fullFileName = fullfile(folder, baseFile);
        if bverbose, fprintf(1, 'Now reading %s\n', fullFileName); end
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
        elseif strfind(baseFile, 'gt')
            s.gtfile = fullFileName;
        elseif strfind(baseFile, 'mask')
            s.maskfile = fullFileName;
        end
    end
    
    %% Normalize intensities
    
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
        
        [ints3d.t1_s2ints, ints3d.t2_s2ints, ints3d.flair_s2ints] = normalizeWM(getfield(load_nii(s.t2_s1file),'img'),...
            getfield(load_nii(s.t2_s2file),'img'),...
            getfield(load_nii(s.flair_s2file),'img'), mask);
    end
    %% distances + blur
    %     dists.t1 = smooth3(ints3d.t1_s2ints - ints3d.t1_s1ints, 'gaussian');
    dists.t1 = ints3d.t1_s2ints - ints3d.t1_s1ints;
    dists.t2 = smooth3(ints3d.t2_s2ints - ints3d.t2_s1ints, 'gaussian');
    dists.flair = ints3d.flair_s2ints - ints3d.flair_s1ints;
    %     dists.flair = smooth3(ints3d.flair_s2ints - ints3d.flair_s1ints, 'gaussian');
    % size(dists.flair)
    %% Voxel selection based on T2-weighted subtraction values
    % dist.t2_dist = ints.t2_s2ints - ints.t2_s1ints;
    
    voxel_selection_mask = zeros(size(dists.t2));
    % size(voxel_selection_mask)
    % numel(smoothedt2_dist)
    % numel(smoothedt2_dist > std(smoothedt2_dist))
    % % std(dists.t2)
    % max(dists.t2(:))
    % numel(voxel_selection_mask(smoothedt2_dist > std(smoothedt2_dist)))
    % std(dists.t2(:))
    %     voxel_selection_mask(dists.t2 > std(dists.t2(mask))) = dists.t2(dists.t2 > std(dists.t2(mask)));
    voxel_selection_mask(dists.t2 > std(dists.t2(mask))) =1;
    
    candidates = voxel_selection_mask;
    % size(candidates)
    % ViewerGUI(candidates)
% %     gt(voxel_selection_mask == 0) = 0;
    %% LRM with Given Parameters
    coefficients = [-9.1008; 0.0021; 0.7388; -0.0540; 0.0001;...
                    -0.2531; 0.6503; -0.0003; 0.5098; -0.8282; 0.0020];
    % size(coefficients)
    DM = @(coef, tdist, t1_s1, t1_s2, t2_s1, t2_s2, flair_s1, flair_s2)...
        coef(1) + (coef(2) * tdist) + (coef(3)*flair_s1) + ...
        (coef(4)*(flair_s2 - flair_s1)) + (coef(5)*tdist*(flair_s2 - flair_s1)) + ...
        (coef(6) * t2_s1) + (coef(7)*(t2_s2 - t2_s1)) + ...
        (coef(8)*tdist*(t2_s2 - t2_s1)) +(coef(9) * t1_s1) + ...
        (coef(10)*(t1_s2 - t1_s1)) + (coef(11)*tdist*(t1_s2 - t1_s1));
    
    pred = DM(coefficients, 0, ints3d.t1_s1ints, ...
        ints3d.t1_s2ints, ints3d.t2_s1ints, ints3d.t2_s2ints, ints3d.flair_s1ints, ints3d.flair_s2ints);
    pred = 1/(1+exp(-pred));
    pred = candidates .* pred;
    %     pred(pred < 0.5) = 0;
    %     pred(pred >= 0.5) = 1;
    %     predictions.% size(test)
    % size(pred)
    % size(voxel_selection_mask)
    %% create predictions
    new = zeros(size(mask));
    % size(new)
    % gaussFilter = gausswin(3);
    % gaussFilter = gaussFilter / sum(gaussFilter);
    % smoothedpred = conv(pred, gaussFilter);
    % size(new(mask))
    % size(smoothedpred)
    % size(mask)
    %     %% save UNCOMMENT LATERAUC
    new(mask) = (pred(mask));
    preds{J} = new;
    
    %     pred(pred < 0) = 0;
    %     pred(pred >= 0) = 1;
    %     temp = make_nii(new);
    %     save_nii(temp, strcat('sweeney_predictions/', 'predictions_', num2str(J),'.nii'));
    %right now this lrm is not returning 0 or 1?
    %% validate
%     gts{J} = gt;
%     [gt, mask] = combine(J, false, false);
    %%
    
        gt = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient', num,'_gt3.nii')),'img');
        gt(gt==2) = 0;
        % plot the ROC curve:
        subplot(1,2,1)
        [X,Y,T,AUCROC,OPTROCPT] = perfcurve(gt(mask>0), new(mask), 1,'xCrit','FPR');%testGroundTruth{J}
        AUCROC
        rocauc(J) = AUCROC
        testAUCROC(s.patient_number) = AUCROC;
        hold on;
        plot(X,Y);
    %     title = strcat
        title('Test Data ROC');
        axis square

        % plot the ROC curve:
        subplot(1,2,2)
        [X,Y,T,AUCPR,OPTROCPT] = perfcurve(gt(mask>0), new(mask), 1,'xCrit','PPV');%testGroundTruth{J
        AUCPR
        prauc(J) = AUCPR
        testAUCPR(s.patient_number) = AUCPR;
        hold on;
        plot(X,Y);
    %     title = strcat
        title('Test Data PR');
        axis square
    
    %     AUC
    % size(new)
    % ViewerGUI(new);
    
    % figure;end
   %% create the volume
%     testpred3d = ind2sub(testVolume(J),testpred);
    dummy_nii = load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_study1_T2Wreg.nii.gz'));
    pred_vol = new;

    %% save predictions as .nii
    temp = make_nii(pred_vol);    
    temp.hdr = dummy_nii.hdr;
    %temp.hdr.dime.bitpix = 256; % signed char
    save_nii(temp, strcat('../MSpatientdata/patient', num, '/patient', num, '_manualDectBaselinetest.nii.gz'));
end
toc


%% Visualize with threshold uncomment later
% N = 5.6;
% new
figure;
thr = 0.2;

seg = new>thr;
% % nnz(seg)
% seg
islice = 26;
subplot(1,2,1);
imagesc(seg(:,:,islice), [0, 1]);
axis image;
% colorbar;
subplot(1,2,2);
imagesc(gt(:,:,islice), [0, 1]);
axis image;
%% validation
% gt(gt<0) = 0;
% validate(new, gt, 20);