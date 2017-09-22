%  loading data into variables for a patient to create structure
patient_number = 1; %change the particular patient's number
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
    elseif strfind(baseFile, 'gt')
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
% s
% s = {t1_s1, t1_s2; t2_s2, t2_s2; flair_s1, flair_s2; mask, gt};
% s.t1_s1 = t1_s1;
%% compute distances
% imagesc(s1t1);
% for sname = {'FLAIR', 'T2'}
%     sname = sname{1};
%     pat.struct.sname;
% end
% dist = zeros(size(t1_s1,1),size(t1_s1,2),size(t1_s1,3));
% dist_t2 = zeros(size(t2_s1,1),size(t2_s1,2),size(t2_s1,3));
% dist_flair = zeros(size(flair_s1,1),size(flair_s1,2),size(flair_s1,3));
% dist = {dist_t1, dist_t2, dist_flair};

% find info about file. need to load before iterating through?
% nii = load_nii(t1_s1);
% A = nii.img;
% size(A) %a x b size images, c slices, d frames. diff b/w frames and slices? want sum over frames/
%53 slices?

% for N = 1:3
%     img1 = s(N, 1);
%     nii = load_nii(img1);
%     img2 = s(N, 2);
%     nii = load_nii(img2);
%     for ii=1:size(img1,1)
%         for jj=1:size(img1,2)
%             for kk=1:size(img2,3)
%                 img2(ii,jj,kk)
%                 %                 distance = (img2(ii,jj,kk) - img2(ii,jj,kk))^2;
%                 %                 dist(ii,jj,kk) = dist(ii,jj,kk) + distance;
%             end
%         end
%     end
% end
% need to initialize as nii load nii first?


% dist = {dist_t1, dist_t2, dist_flair};
% dist
%
%% Generate matrices of intensities

nii = load_nii(s.maskfile);
mask = nii.img;
mask = logical(mask);
nii = load_nii(s.gtfile);
gt = nii.img;
% gt;
fields = fieldnames(s);

for N = 2:numel(fields) %skip patient number, go directly to index 2.
    field = fields{N};
    value = getfield(s, field);
    nii = load_nii(value);
    im = nii.img;
    r= strrep(field,'file','');
    ints.(strcat(r,'ints')) = im(mask);
end
% intensities
%% Find Euclidean distance
dist.t1_dist = (ints.t1_s2ints - ints.t1_s1ints).^2;
dist.t2_dist = (ints.t2_s2ints - ints.t2_s1ints).^2;
dist.flair_dist= (ints.flair_s2ints - ints.flair_s1ints).^2;
dist.total = (dist.t1_dist + dist.t2_dist + dist.flair_dist).^0.5;
%% Visualize
new = zeros(size(mask));
new(mask) = dist.total;
ViewerGUI(new);
%% validate

validate(new, gt, 10);
% %% Validation:sensitivity, specificity, accuracy, positive predictive value, and dice similarity index
% %sensitivity: amount deteceted. specificity: true negative/all negatives
% %accuracies: al true/all ppv = true postiive/all postive dsi = 2 * true
% %positives / each of segmented sets. rock curve, PR curve
% % use parameters seg + gt(mask) ? or just gt.
% % make sure both are boolean matrices
% se = @(auto_seg,manual_seg) nnz(auto_seg & manual_seg)/nnz(manual_seg); %here tp/fn + tp = recall
% sp = @(auto_seg,manual_seg) nnz(~auto_seg & ~manual_seg)/nnz(~manual_seg); %denominator is all negatives detected
% a = @(auto_seg,manual_seg) nnz(or(~auto_seg & ~manual_seg, auto_seg & manual_seg))/numel(manual_seg);
% ppv = @(auto_seg,manual_seg) nnz(auto_seg & manual_seg)/nnz(auto_seg); %here tp/fp + tp = precision
% dice = @(auto_seg,manual_seg) 2*nnz(auto_seg & manual_seg)/(nnz(auto_seg) + nnz(manual_seg));
%% Visualize with threshold
N = 2;
% new
thr = 150;
seg = new>thr;
% nnz(seg)n
% seg
% figure;
islice = 25;
subplot(1,2,1);
imagesc(seg(:,:,islice), [0, 1]);
axis image;
% colorbar;
subplot(1,2,2);
imagesc(gt(:,:,islice), [0, 1]);
axis image;
% end
% gt

%% Select slice and visualize that slice
islice = 37;
imagesc(new(:,:,islice), [min(dist.total), max(dist.total)]);
axis image;
colorbar;
%% Select slice and visualize that slice + side by side with ground truth
islice = 26;
subplot(1,2,1);
imagesc(new(:,:,islice), [min(dist.total), max(dist.total)]);
axis image;
% colorbar;
subplot(1,2,2);
imagesc(gt(:,:,islice), [0, 1]);
axis image;
%% Select slice and visualize that slice with gt overlap
islice = 31;
clf;
imagesc(new(:,:,islice), [min(dist.total), max(dist.total)]);
axis image;
colorbar;
hold on;
contour(gt(:,:,islice), 1, 'r');