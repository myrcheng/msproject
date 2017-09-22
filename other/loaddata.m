function [data, gt] = loaddata(patient_number)
%% loading data into variables for a patient to create structure
% patient_number = 1; %change the particular patient's number
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
end

%% Generate intensities

nii = load_nii(s.maskfile);
mask = nii.img;
mask = logical(mask);
nii = load_nii(s.gtfile);
gt = nii.img;
fields = fieldnames(s);

for N = 2:numel(fields) %skip patient number, go directly to index 2.
    field = fields{N};
    value = getfield(s, field);
    nii = load_nii(value);
    im = nii.img;
    r= strrep(field,'file','');
    temp = im(mask);
%     a = std(temp(:))
%        mean(temp);

    centered = (temp - mean(temp));
    sum = 0;
    for i = 1:numel(centered)
        sum = sum + (centered(i, 1)^2);
    end
    stddev = sqrt(double(sum)/length(centered));
    ints3d.(strcat(r,'ints')) = zeros(size(mask));
    ints3d.(strcat(r,'ints'))(mask) = centered/stddev; %1d matrix
end
%% distances + blur
dists.t1 = smooth3(ints3d.t1_s2ints - ints3d.t1_s1ints, 'gaussian', 5);
dists.t2 = smooth3(ints3d.t2_s2ints - ints3d.t2_s1ints, 'gaussian', 5);
dists.flair = smooth3(ints3d.flair_s2ints - ints3d.flair_s1ints, 'gaussian', 5);
dfields = fieldnames(dists);
% for i = 1:numel(dfields)
%     field = dfields{i};
%     value = getfield(dists, field)(mask);
%     dists.{field = value;
% end
dists.t1 = dists.t1(mask);
dists.t2 = dists.t2(mask);
dists.flair = dists.flair(mask);

%% Voxel selection based on T2-weighted subtraction values

voxel_selection_mask = zeros(size(ints3d.flair_s2ints(mask)));
size(voxel_selection_mask)
size(dists.t2(dists.t2 > std(dists.t2(:))))
voxel_selection_mask(dists.t2 > std(dists.t2(:))) = dists.t2(dists.t2 > std(dists.t2(:)));
inds = find(mask);
gt(inds(dists.t2 <= std(dists.t2(:)))) = 0;
se = nnz(voxel_selection_mask & (gt(mask) > 0))/nnz(gt(mask)>0)
% candidates = voxel_selection_mask;
% gt
%% Apply voxel selection mask
dists.t1(voxel_selection_mask == 0) = [];
dists.t2(voxel_selection_mask == 0) = [];
dists.flair(voxel_selection_mask == 0) = [];
ints3d.flair_s2ints = ints3d.flair_s2ints(mask);
ints3d.t1_s2ints = ints3d.t1_s2ints(mask);
ints3d.t2_s2ints = ints3d.t2_s2ints(mask);
ints3d.flair_s2ints(voxel_selection_mask == 0) = [];
ints3d.t1_s2ints(voxel_selection_mask == 0) = [];
ints3d.t2_s2ints(voxel_selection_mask == 0) = [];
gt = gt(mask);
gt(voxel_selection_mask == 0) = [];
% size(tempdata)
% size(data)
% for bridge = each(fieldnames(signal))
%    signal.(bridge) = rand(10);
% end
%% Load data into prediction
deltaT = ones(size(dists.flair));
deltaT = 365 * deltaT;
data = [ints3d.flair_s2ints, ints3d.t2_s2ints, ints3d.t1_s2ints, dists.flair, dists.t2, dists.t1, deltaT, (deltaT .* dists.flair), (deltaT .* dists.t2), (deltaT .* dists.t1)];
% data = [ints3d.flair_s2ints(mask), ints3d.t2_s2ints(mask), ints3d.t1_s2ints(mask), dists.flair(mask), dists.t2(mask), dists.t1(mask), deltaT, (deltaT .* dists.flair(mask)), (deltaT .* dists.t2(mask)), (deltaT .* dists.t1(mask))];
% tempdata = data;
size(data)
% data(data < 0) = 0;
end