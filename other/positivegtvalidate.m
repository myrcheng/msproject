gts = cell(20, 1);
auc = zeros(20, 1);
for i = 1:20
    %     init/load nii and generate ground truth
    num = num2str(i);
    fileName = strcat('../MSproject/sweeney_predictions/predictions_', num, '.nii');
    fprintf(1, 'Now reading %s\n', fileName);
    [gt, mask] = combine(i, false);
    nii = load_nii(fileName);
    new = nii.img;
    %     new = logical(mask);
    %     new(mask) = mask;
    %     mask = logical(mask);
    %     %% Validation:sensitivity, specificity, accuracy, positive predictive value, and dice similarity index
    %     %sensitivity: amount deteceted. specificity: true negative/all negatives
    %     %accuracies: al true/all ppv = true postiive/all postive dsi = 2 * true
    %     %positives / each of segmented sets. rock curve, PR curve
    %     % use parameters seg + gt(mask) ? or just gt.
    %     % make sure both are boolean matrices
    %     se = @(auto_seg,manual_seg) nnz(auto_seg & manual_seg)/nnz(manual_seg); %here tp/fn + tp = recall
    %     sp = @(auto_seg,manual_seg) nnz(~auto_seg & ~manual_seg)/nnz(~manual_seg); %denominator is all negatives detected
    %     a = @(auto_seg,manual_seg) nnz(or(~auto_seg & ~manual_seg, auto_seg & manual_seg))/numel(manual_seg);
    %     ppv = @(auto_seg,manual_seg) nnz(auto_seg & manual_seg)/nnz(auto_seg); %here tp/fp + tp = precision
    %     dice = @(auto_seg,manual_seg) 2*nnz(auto_seg & manual_seg)/(nnz(auto_seg) + nnz(manual_seg));
    %     %% Generate table of data
    %     pts = 10;
    %     %only care about masked part...try more though just to graph it.
    %     sem = zeros(pts, 1);
    %     spm = zeros(pts, 1);
    %     am = zeros(pts, 1);
    %     ppvm = zeros(pts, 1);
    %     dicem = zeros(pts, 1);
    %
    %     % sem, spm, am, ppvm, dicem = deal(zeros(10,1));
    %     for N = 1:pts
    %         thr = N;
    %         seg = new>thr;
    %         gttest = (gt == 1);
    %         sem(N) = se(seg, gttest);
    %         spm(N) = sp(seg, gttest);
    %         am(N) = a(seg, gttest);
    %         ppvm(N) = ppv(seg, gttest);
    %         dicem(N) = dice(seg, gttest);
    %
    %     end
    %% Plots of validation measurements UNCOMMENT
    %     figure;
    %     plot(sem, '-or');
    %     hold on;
    %     plot(spm, '-om');
    %     hold on;
    %     plot(am, '-oc');
    %     hold on;
    %     plot(ppvm, '-og');
    %     hold on;
    %     plot(dicem, '-ob');
    %     hold on;
    % plot
    %% Plots of curves - PR and AUC UNCOMMENT
    % precision = ppv, recall = sv
    % 1 - sp = fpr = fp / n
    %     figure;
    %     plot(sem, ppvm);
    %     title('PR Curve');
    %     xlabel('Recall');
    %     ylabel('Precision');
    %     figure;
    %     plot(1-spm, sem, '-or');
    %     % 1-spm
    %     title('ROC Curve');
    %     xlabel('False Positive Rate (1-specificity)');
    %     ylabel('Sensitivity');
    %     spm = [spm; 0; 1];
    %     sem = [sem; 1; 0];
    % auc = trapz(1-spm, sem)
    %% perfcurve function
    gts{i} = gt;
    if any(gt(:) == 0 | gt(:) == 2) && any(gt(:) == 1)
        [X,Y,T,AUC] = perfcurve(gt(mask) == 1, new(mask), true);
        %new(mask) just shows whether it is part of the brain mask or not
        auc(i) = AUC;
        AUC
    end
end