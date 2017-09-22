function [feat,featureGroupID, featureGroupList, featureList] = computeFeats(img, roi)
% roi = true(21);
roi = ones(21);
if nargin<2, roi = ones(size(img)); end
GrayLimits = [-6,6];

featureGroupList = {'Histogram';...
    'Intensity Difference';...
    'Elliptical Fit';...
    'GLCM';...
    'GLRLM'};
i = 0;
igroup = 0;
% 1. histogram group:
igroup = igroup + 1;
i = i+1;
feat(i) = mean(img(:));

% feat(i) = mean(img(roi));
featureList{i} = 'Mean';
featureGroupID(i) = igroup;

i = i+1;
feat(i) = std(img(:));
% feat(i) = std(img(roi));
featureList{i} = 'StandardDeviation';
featureGroupID(i) = igroup;

i = i+1;
feat(i) = skewness(img(:));
featureList{i} = 'Skewness';
featureGroupID(i) = igroup;

i = i+1;
feat(i) = kurtosis(img(:));
featureList{i} = 'Kurtosis';
featureGroupID(i) = igroup;

% 4. GLCM:
igroup = igroup + 1;

tmp = double(roi>0);
tmp(roi==0) = min(GrayLimits)-1;
glcm2 = graycomatrix(img.*tmp,'Offset',[2 0;0 2],'NumLevels',5,'GrayLimits',GrayLimits);
glcm = mean(glcm2,3);
s = GLCM_Features1(glcm);

names = {'Autocorrelation','Contrast','Correlation','ClusterProminence',...
    'ClusterShade', 'Dissimilarity', 'Energy', 'Entropy',...
    'Homogeneity', 'Variance', 'SumAverage',...
    'MaximumProbability', 'SumVariance', 'SumEntropy',...
    'DifferenceVariance','DifferenceEntropy'};
values = {s.autoc, s.contr, s.corrm, s.cprom,...
    s.cshad, s.dissi, s.energ, s.entro,...
    s.homom, s.sosvh, s.savgh,...
    s.maxpr, s.svarh, s.senth, s.dvarh, s.denth};
for ii = 1:numel(names)
    i = i+1;
    feat(i) = values{ii};
    featureList{i} = names{ii};
    featureGroupID(i) = igroup;
end

% 5. GLRLM
igroup = igroup + 1;

tmp = double(roi>0);
tmp(roi==0) = min(GrayLimits)-1;
glrlms = grayrlmatrix(img.*tmp, 'GrayLimits',GrayLimits);
S = grayrlprops(glrlms);
s = struct;
snames = fieldnames(S);
for sfield = snames(1:numel(snames)-1)'
    s.(sfield{1}) = mean(S.(sfield{1}));
end
names = {'SRE', 'LRE', 'GLN', 'RLN', 'RP', 'LGRE',...
    'HGRE', 'SGLGE', 'SRHGE', 'LRLGE', 'LRHGE'};
values = {s.SRE, s.LRE, s.GLN, s.RLN, s.RP, s.LGRE,...
    s.HGRE, s.SGLGE, s.SRHGE, s.LRLGE, s.LRHGE};
for ii = 1:numel(names)
    i = i+1;
    feat(i) = values{ii};
    featureList{i} = names{ii};
    featureGroupID(i) = igroup;
end

feat = feat';
featureGroupID = featureGroupID';
featureList = featureList';

%% %%%%%%%%%%%%%%%%%%%%%%%% helpers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = GLCM_Features1(glcmin,pairs)
% 
% GLCM_Features1 helps to calculate the features from the different GLCMs
% that are input to the function. The GLCMs are stored in a i x j x n
% matrix, where n is the number of GLCMs calculated usually due to the
% different orientation and displacements used in the algorithm. Usually
% the values i and j are equal to 'NumLevels' parameter of the GLCM
% computing function graycomatrix(). Note that matlab quantization values
% belong to the set {1,..., NumLevels} and not from {0,...,(NumLevels-1)}
% as provided in some references
% http://www.mathworks.com/access/helpdesk/help/toolbox/images/graycomatrix
% .html
% 
% Although there is a function graycoprops() in Matlab Image Processing
% Toolbox that computes four parameters Contrast, Correlation, Energy,
% and Homogeneity. The paper by Haralick suggests a few more parameters
% that are also computed here. The code is not fully vectorized and hence
% is not an efficient implementation but it is easy to add new features
% based on the GLCM using this code. Takes care of 3 dimensional glcms
% (multiple glcms in a single 3D array)
% 
% If you find that the values obtained are different from what you expect 
% or if you think there is a different formula that needs to be used 
% from the ones used in this code please let me know. 
% A few questions which I have are listed in the link 
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/239608
%
% I plan to submit a vectorized version of the code later and provide 
% updates based on replies to the above link and this initial code. 
%
% Features computed 
% Autocorrelation: [2]                      (out.autoc)
% Contrast: matlab/[1,2]                    (out.contr)
% Correlation: matlab                       (out.corrm)
% Correlation: [1,2]                        (out.corrp)
% Cluster Prominence: [2]                   (out.cprom)
% Cluster Shade: [2]                        (out.cshad)
% Dissimilarity: [2]                        (out.dissi)
% Energy: matlab / [1,2]                    (out.energ)
% Entropy: [2]                              (out.entro)
% Homogeneity: matlab                       (out.homom)
% Homogeneity: [2]                          (out.homop)
% Maximum probability: [2]                  (out.maxpr)
% Sum of sqaures: Variance [1]              (out.sosvh)
% Sum average [1]                           (out.savgh)
% Sum variance [1]                          (out.svarh)
% Sum entropy [1]                           (out.senth)
% Difference variance [1]                   (out.dvarh)
% Difference entropy [1]                    (out.denth)
% Information measure of correlation1 [1]   (out.inf1h)
% Informaiton measure of correlation2 [1]   (out.inf2h)
% Inverse difference (INV) is homom [3]     (out.homom)
% Inverse difference normalized (INN) [3]   (out.indnc) 
% Inverse difference moment normalized [3]  (out.idmnc)
%
% The maximal correlation coefficient was not calculated due to
% computational instability 
% http://murphylab.web.cmu.edu/publications/boland/boland_node26.html
%
% Formulae from MATLAB site (some look different from
% the paper by Haralick but are equivalent and give same results)
% Example formulae: 
% Contrast = sum_i(sum_j(  (i-j)^2 * p(i,j) ) ) (same in matlab/paper)
% Correlation = sum_i( sum_j( (i - u_i)(j - u_j)p(i,j)/(s_i.s_j) ) ) (m)
% Correlation = sum_i( sum_j( ((ij)p(i,j) - u_x.u_y) / (s_x.s_y) ) ) (p[2])
% Energy = sum_i( sum_j( p(i,j)^2 ) )           (same in matlab/paper)
% Homogeneity = sum_i( sum_j( p(i,j) / (1 + |i-j|) ) ) (as in matlab)
% Homogeneity = sum_i( sum_j( p(i,j) / (1 + (i-j)^2) ) ) (as in paper)
% 
% Where:
% u_i = u_x = sum_i( sum_j( i.p(i,j) ) ) (in paper [2])
% u_j = u_y = sum_i( sum_j( j.p(i,j) ) ) (in paper [2])
% s_i = s_x = sum_i( sum_j( (i - u_x)^2.p(i,j) ) ) (in paper [2])
% s_j = s_y = sum_i( sum_j( (j - u_y)^2.p(i,j) ) ) (in paper [2])
%
% 
% Normalize the glcm:
% Compute the sum of all the values in each glcm in the array and divide 
% each element by it sum
%
% Haralick uses 'Symmetric' = true in computing the glcm
% There is no Symmetric flag in the Matlab version I use hence
% I add the diagonally opposite pairs to obtain the Haralick glcm
% Here it is assumed that the diagonally opposite orientations are paired
% one after the other in the matrix
% If the above assumption is true with respect to the input glcm then
% setting the flag 'pairs' to 1 will compute the final glcms that would result 
% by setting 'Symmetric' to true. If your glcm is computed using the
% Matlab version with 'Symmetric' flag you can set the flag 'pairs' to 0
%
% References:
% 1. R. M. Haralick, K. Shanmugam, and I. Dinstein, Textural Features of
% Image Classification, IEEE Transactions on Systems, Man and Cybernetics,
% vol. SMC-3, no. 6, Nov. 1973
% 2. L. Soh and C. Tsatsoulis, Texture Analysis of SAR Sea Ice Imagery
% Using Gray Level Co-Occurrence Matrices, IEEE Transactions on Geoscience
% and Remote Sensing, vol. 37, no. 2, March 1999.
% 3. D A. Clausi, An analysis of co-occurrence texture statistics as a
% function of grey level quantization, Can. J. Remote Sensing, vol. 28, no.
% 1, pp. 45-62, 2002
% 4. http://murphylab.web.cmu.edu/publications/boland/boland_node26.html
%
%
% Example:
%
% Usage is similar to graycoprops() but needs extra parameter 'pairs' apart
% from the GLCM as input
% I = imread('circuit.tif');
% GLCM2 = graycomatrix(I,'Offset',[2 0;0 2]);
% stats = GLCM_features1(GLCM2,0)
% The output is a structure containing all the parameters for the different
% GLCMs
%
% [Avinash Uppuluri: avinash_uv@yahoo.com: Last modified: 11/20/08]

% If 'pairs' not entered: set pairs to 0 
if ((nargin > 2) || (nargin == 0))
   error('Too many or too few input arguments. Enter GLCM and pairs.');
elseif ( (nargin == 2) ) 
    if ((size(glcmin,1) <= 1) || (size(glcmin,2) <= 1))
       error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcmin,1) ~= size(glcmin,2) )
        error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end    
elseif (nargin == 1) % only GLCM is entered
    pairs = 0; % default is numbers and input 1 for percentage
    if ((size(glcmin,1) <= 1) || (size(glcmin,2) <= 1))
       error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcmin,1) ~= size(glcmin,2) )
       error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end    
end


format long e
if (pairs == 1)
    newn = 1;
    for nglcm = 1:2:size(glcmin,3)
        glcm(:,:,newn)  = glcmin(:,:,nglcm) + glcmin(:,:,nglcm+1);
        newn = newn + 1;
    end
elseif (pairs == 0)
    glcm = glcmin;
end

size_glcm_1 = size(glcm,1);
size_glcm_2 = size(glcm,2);
size_glcm_3 = size(glcm,3);

% checked 
out.autoc = zeros(1,size_glcm_3); % Autocorrelation: [2] 
out.contr = zeros(1,size_glcm_3); % Contrast: matlab/[1,2]
out.corrm = zeros(1,size_glcm_3); % Correlation: matlab
out.corrp = zeros(1,size_glcm_3); % Correlation: [1,2]
out.cprom = zeros(1,size_glcm_3); % Cluster Prominence: [2]
out.cshad = zeros(1,size_glcm_3); % Cluster Shade: [2]
out.dissi = zeros(1,size_glcm_3); % Dissimilarity: [2]
out.energ = zeros(1,size_glcm_3); % Energy: matlab / [1,2]
out.entro = zeros(1,size_glcm_3); % Entropy: [2]
out.homom = zeros(1,size_glcm_3); % Homogeneity: matlab
out.homop = zeros(1,size_glcm_3); % Homogeneity: [2]
out.maxpr = zeros(1,size_glcm_3); % Maximum probability: [2]

out.sosvh = zeros(1,size_glcm_3); % Sum of sqaures: Variance [1]
out.savgh = zeros(1,size_glcm_3); % Sum average [1]
out.svarh = zeros(1,size_glcm_3); % Sum variance [1]
out.senth = zeros(1,size_glcm_3); % Sum entropy [1]
out.dvarh = zeros(1,size_glcm_3); % Difference variance [4]
%out.dvarh2 = zeros(1,size_glcm_3); % Difference variance [1]
out.denth = zeros(1,size_glcm_3); % Difference entropy [1]
out.inf1h = zeros(1,size_glcm_3); % Information measure of correlation1 [1]
out.inf2h = zeros(1,size_glcm_3); % Informaiton measure of correlation2 [1]
%out.mxcch = zeros(1,size_glcm_3);% maximal correlation coefficient [1]
%out.invdc = zeros(1,size_glcm_3);% Inverse difference (INV) is homom [3]
out.indnc = zeros(1,size_glcm_3); % Inverse difference normalized (INN) [3]
out.idmnc = zeros(1,size_glcm_3); % Inverse difference moment normalized [3]

% correlation with alternate definition of u and s
%out.corrm2 = zeros(1,size_glcm_3); % Correlation: matlab
%out.corrp2 = zeros(1,size_glcm_3); % Correlation: [1,2]

glcm_sum  = zeros(size_glcm_3,1);
glcm_mean = zeros(size_glcm_3,1);
glcm_var  = zeros(size_glcm_3,1);

% http://www.fp.ucalgary.ca/mhallbey/glcm_mean.htm confuses the range of 
% i and j used in calculating the means and standard deviations.
% As of now I am not sure if the range of i and j should be [1:Ng] or
% [0:Ng-1]. I am working on obtaining the values of mean and std that get
% the values of correlation that are provided by matlab.
u_x = zeros(size_glcm_3,1);
u_y = zeros(size_glcm_3,1);
s_x = zeros(size_glcm_3,1);
s_y = zeros(size_glcm_3,1);

% % alternate values of u and s
% u_x2 = zeros(size_glcm_3,1);
% u_y2 = zeros(size_glcm_3,1);
% s_x2 = zeros(size_glcm_3,1);
% s_y2 = zeros(size_glcm_3,1);

% checked p_x p_y p_xplusy p_xminusy
p_x = zeros(size_glcm_1,size_glcm_3); % Ng x #glcms[1]  
p_y = zeros(size_glcm_2,size_glcm_3); % Ng x #glcms[1]
p_xplusy = zeros((size_glcm_1*2 - 1),size_glcm_3); %[1]
p_xminusy = zeros((size_glcm_1),size_glcm_3); %[1]
% checked hxy hxy1 hxy2 hx hy
hxy  = zeros(size_glcm_3,1);
hxy1 = zeros(size_glcm_3,1);
hx   = zeros(size_glcm_3,1);
hy   = zeros(size_glcm_3,1);
hxy2 = zeros(size_glcm_3,1);

%Q    = zeros(size(glcm));

for k = 1:size_glcm_3 % number glcms

    glcm_sum(k) = sum(sum(glcm(:,:,k)));
    glcm(:,:,k) = glcm(:,:,k)./glcm_sum(k); % Normalize each glcm
    glcm_mean(k) = mean2(glcm(:,:,k)); % compute mean after norm
    glcm_var(k)  = (std2(glcm(:,:,k)))^2;
    
    for i = 1:size_glcm_1

        for j = 1:size_glcm_2

            out.contr(k) = out.contr(k) + (abs(i - j))^2.*glcm(i,j,k);
            out.dissi(k) = out.dissi(k) + (abs(i - j)*glcm(i,j,k));
            out.energ(k) = out.energ(k) + (glcm(i,j,k).^2);
            out.entro(k) = out.entro(k) - (glcm(i,j,k)*log(glcm(i,j,k) + eps));
            out.homom(k) = out.homom(k) + (glcm(i,j,k)/( 1 + abs(i-j) ));
            out.homop(k) = out.homop(k) + (glcm(i,j,k)/( 1 + (i - j)^2));
            % [1] explains sum of squares variance with a mean value;
            % the exact definition for mean has not been provided in 
            % the reference: I use the mean of the entire normalized glcm 
            out.sosvh(k) = out.sosvh(k) + glcm(i,j,k)*((i - glcm_mean(k))^2);
            
            %out.invdc(k) = out.homom(k);
            out.indnc(k) = out.indnc(k) + (glcm(i,j,k)/( 1 + (abs(i-j)/size_glcm_1) ));
            out.idmnc(k) = out.idmnc(k) + (glcm(i,j,k)/( 1 + ((i - j)/size_glcm_1)^2));
            u_x(k)          = u_x(k) + (i)*glcm(i,j,k); % changed 10/26/08
            u_y(k)          = u_y(k) + (j)*glcm(i,j,k); % changed 10/26/08
            % code requires that Nx = Ny 
            % the values of the grey levels range from 1 to (Ng) 
        end
        
    end
    out.maxpr(k) = max(max(glcm(:,:,k)));
end
% glcms have been normalized:
% The contrast has been computed for each glcm in the 3D matrix
% (tested) gives similar results to the matlab function

for k = 1:size_glcm_3
    
    for i = 1:size_glcm_1
        
        for j = 1:size_glcm_2
            p_x(i,k) = p_x(i,k) + glcm(i,j,k); 
            p_y(i,k) = p_y(i,k) + glcm(j,i,k); % taking i for j and j for i
            if (ismember((i + j),[2:2*size_glcm_1])) 
                p_xplusy((i+j)-1,k) = p_xplusy((i+j)-1,k) + glcm(i,j,k);
            end
            if (ismember(abs(i-j),[0:(size_glcm_1-1)])) 
                p_xminusy((abs(i-j))+1,k) = p_xminusy((abs(i-j))+1,k) +...
                    glcm(i,j,k);
            end
        end
    end
    
%     % consider u_x and u_y and s_x and s_y as means and standard deviations
%     % of p_x and p_y
%     u_x2(k) = mean(p_x(:,k));
%     u_y2(k) = mean(p_y(:,k));
%     s_x2(k) = std(p_x(:,k));
%     s_y2(k) = std(p_y(:,k));
    
end

% marginal probabilities are now available [1]
% p_xminusy has +1 in index for matlab (no 0 index)
% computing sum average, sum variance and sum entropy:
for k = 1:(size_glcm_3)
    
    for i = 1:(2*(size_glcm_1)-1)
        out.savgh(k) = out.savgh(k) + (i+1)*p_xplusy(i,k);
        % the summation for savgh is for i from 2 to 2*Ng hence (i+1)
        out.senth(k) = out.senth(k) - (p_xplusy(i,k)*log(p_xplusy(i,k) + eps));
    end

end
% compute sum variance with the help of sum entropy
for k = 1:(size_glcm_3)
    
    for i = 1:(2*(size_glcm_1)-1)
        out.svarh(k) = out.svarh(k) + (((i+1) - out.senth(k))^2)*p_xplusy(i,k);
        % the summation for savgh is for i from 2 to 2*Ng hence (i+1)
    end

end
% compute difference variance, difference entropy, 
for k = 1:size_glcm_3
% out.dvarh2(k) = var(p_xminusy(:,k));
% but using the formula in 
% http://murphylab.web.cmu.edu/publications/boland/boland_node26.html
% we have for dvarh
    for i = 0:(size_glcm_1-1)
        out.denth(k) = out.denth(k) - (p_xminusy(i+1,k)*log(p_xminusy(i+1,k) + eps));
        out.dvarh(k) = out.dvarh(k) + (i^2)*p_xminusy(i+1,k);
    end
end

% compute information measure of correlation(1,2) [1]
for k = 1:size_glcm_3
    hxy(k) = out.entro(k);
    for i = 1:size_glcm_1
        
        for j = 1:size_glcm_2
            hxy1(k) = hxy1(k) - (glcm(i,j,k)*log(p_x(i,k)*p_y(j,k) + eps));
            hxy2(k) = hxy2(k) - (p_x(i,k)*p_y(j,k)*log(p_x(i,k)*p_y(j,k) + eps));
%             for Qind = 1:(size_glcm_1)
%                 Q(i,j,k) = Q(i,j,k) +...
%                     ( glcm(i,Qind,k)*glcm(j,Qind,k) / (p_x(i,k)*p_y(Qind,k)) ); 
%             end
        end
        hx(k) = hx(k) - (p_x(i,k)*log(p_x(i,k) + eps));
        hy(k) = hy(k) - (p_y(i,k)*log(p_y(i,k) + eps));
    end
    out.inf1h(k) = ( hxy(k) - hxy1(k) ) / ( max([hx(k),hy(k)]) );
    out.inf2h(k) = ( 1 - exp( -2*( hxy2(k) - hxy(k) ) ) )^0.5;
%     eig_Q(k,:)   = eig(Q(:,:,k));
%     sort_eig(k,:)= sort(eig_Q(k,:),'descend');
%     out.mxcch(k) = sort_eig(k,2)^0.5;
% The maximal correlation coefficient was not calculated due to
% computational instability 
% http://murphylab.web.cmu.edu/publications/boland/boland_node26.html
end

corm = zeros(size_glcm_3,1);
corp = zeros(size_glcm_3,1);
% using http://www.fp.ucalgary.ca/mhallbey/glcm_variance.htm for s_x s_y
for k = 1:size_glcm_3
    for i = 1:size_glcm_1
        for j = 1:size_glcm_2
            s_x(k)  = s_x(k)  + (((i) - u_x(k))^2)*glcm(i,j,k);
            s_y(k)  = s_y(k)  + (((j) - u_y(k))^2)*glcm(i,j,k);
            corp(k) = corp(k) + ((i)*(j)*glcm(i,j,k));
            corm(k) = corm(k) + (((i) - u_x(k))*((j) - u_y(k))*glcm(i,j,k));
            out.cprom(k) = out.cprom(k) + (((i + j - u_x(k) - u_y(k))^4)*...
                glcm(i,j,k));
            out.cshad(k) = out.cshad(k) + (((i + j - u_x(k) - u_y(k))^3)*...
                glcm(i,j,k));
        end
    end
    % using http://www.fp.ucalgary.ca/mhallbey/glcm_variance.htm for s_x
    % s_y : This solves the difference in value of correlation and might be
    % the right value of standard deviations required 
    % According to this website there is a typo in [2] which provides
    % values of variance instead of the standard deviation hence a square
    % root is required as done below:
    s_x(k) = s_x(k) ^ 0.5;
    s_y(k) = s_y(k) ^ 0.5;
    out.autoc(k) = corp(k);
    out.corrp(k) = (corp(k) - u_x(k)*u_y(k))/(s_x(k)*s_y(k));
    out.corrm(k) = corm(k) / (s_x(k)*s_y(k));
%     % alternate values of u and s
%     out.corrp2(k) = (corp(k) - u_x2(k)*u_y2(k))/(s_x2(k)*s_y2(k));
%     out.corrm2(k) = corm(k) / (s_x2(k)*s_y2(k));
end

function [GLRLMS,SI]= grayrlmatrix(varargin)
%  Description
%  -------------------------------------------
%   Computes the graylevel run length (GLRL) matrix used for textural
%   analysis of an image using zigzag scan method.The method includes four
%   basic steps
%       Step 1 determine direction
%       Step 2 zigzag scan
%       Step 3 obtain new sequences
%       Step 4 calculate run-length matrix
%   -----------------------------------------
%   GLRLMS = GRAYRLMATRIX(I,PARAM1,VALUE1,PARAM2,VALUE2,...) returns one or more
%   gray-level run-length matrices, depending on the values of the optional
%   parameter/value pairs. Parameter names can be abbreviated, and case does
%   not matter.
%  ------------------------------------------
%   Parameters include:
%  ------------------------------------------
%   'Offset'         A p-by-1 vector of offsets specifying the scanning direction.
%
%
%                    Angle     OFFSET
%                    -----     ------
%                    0          1
%                    45         2
%                    90         3
%                    135        4
%
%                    OFFSET must be integers from {1 2 3 4}.
%
%                    Default: [1 2 3 4]
%
%   'NumLevels'      An integer specifying the number of gray levels to use when
%                    scaling the grayscale values in I. For example, if
%                    'NumLevels' is 8, GRAYRLMATRIX scales the values in I so
%                    they are integers between 1 and 8.  The number of gray levels
%                    determines the size of the gray-level run-length matrix
%
%
%                    'NumLevels' must be an integer. 'NumLevels' must be 2 if I
%                    is logical.
%
%                    Default: 8 for numeric
%                             2 for logical
%
%   'GrayLimits'     A two-element vector, [LOW HIGH], that specifies how the
%                    grayscale values in I are linearly scaled into gray
%                    levels. Grayscale values less than or equal to LOW are
%                    scaled to 1. Grayscale values greater than or equal to
%                    HIGH are scaled to HIGH.  If 'GrayLimits' is set to [],
%                    GRAYRLMATRIX uses the minimum and maximum grayscale values
%                    in I as limits, [min(I(:)) max(I(:))].
%
%                    Default: the LOW and HIGH values specified by the
%                    class, e.g., [LOW HIGH] is [0 1] if I is double and
%                    [-32768 32767] if I is int16.
%
%  ------------------------------------------
%  Example
%  ------------------------------------------
% I =[1     1     1     2     2
%      3     4     2     2     3
%      4     4     4     4     4
%      5     5     3     3     3
%      1     1     3     4     5]
% [GLRLMS,SI] = grayrlmatrix(I,'NumLevels',5,'G',[])
% I =
%      1     1     1     2     2
%      3     4     2     2     3
%      4     4     4     4     4
%      5     5     3     3     3
%      1     1     3     4     5
% GLRLMS(:,:,1) =
%      0     1     1     0     0
%      0     2     0     0     0
%      3     0     1     0     0
%      2     0     0     0     1
%      1     1     0     0     0
% GLRLMS(:,:,2) =
%      5     0     0     0     0
%      0     2     0     0     0
%      4     1     0     0     0
%      5     1     0     0     0
%      3     0     0     0     0
% GLRLMS(:,:,3) =
%      5     0     0     0     0
%      2     1     0     0     0
%      4     1     0     0     0
%      5     1     0     0     0
%      3     0     0     0     0
% GLRLMS(:,:,4) =
%      5     0     0     0     0
%      4     0     0     0     0
%      6     0     0     0     0
%      5     1     0     0     0
%      3     0     0     0     0
% SI =
%      1     1     1     2     2
%      3     4     2     2     3
%      4     4     4     4     4
%      5     5     3     3     3
%      1     1     3     4     5
% -------------------------------------------
% See also zigzag rle_0 rle_45
% -------------------------------------------
% Author:
% -------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,100076
% -------------------------------------------
% History:
% -------------------------------------------
% Creation: beta  Date: 01/10/2007
% Revision: 1.0   Date: 14/11/2007
% -------------------------------------------
% Bug Fixed:
% -------------------------------------------
% 1.Issue wrong results for nonsquare matrix,now output cells instead of
%   multi-dim arrays
% 2.Add support for inputs checking inspired by MATLAB style
% 
[I, Offset, NL, GL] = ParseInputs(varargin{:});

% Scale I so that it contains integers between 1 and NL.
if GL(2) == GL(1)
    SI = ones(size(I));
else
    slope = (NL-1) / (GL(2) - GL(1));
    intercept = 1 - (slope*(GL(1)));
    SI = round(imlincomb(slope,I,intercept,'double'));
end

% Clip values if user had a value that is outside of the range, e.g., double
% image = [0 .5 2;0 1 1]; 2 is outside of [0,1]. The order of the following
% lines matters in the event that NL = 0.
SI(SI > NL) = NL;
SI(SI < 1) = 1;
% total numbers of direction
numOffsets = size(Offset,1);

if NL ~= 0
    % make direction matrix for all given directions
    for k = 1 : numOffsets
        GLRLMS{k} = computeGLRLM(SI,Offset(k),NL);
    end
else
    GLRLMS = [];
end

% --------------------------------------------------------------------
function oneGLRLM = computeGLRLM(si,offset,nl)
% For given direction, compute the run length matrix
switch offset
    case 1
        % 0 degree
        oneGLRLM = rle_0(si,nl);
    case 2
        % 45 degree
        seq = zigzag(si);
        oneGLRLM  = rle_45(seq,nl);
    case 3
        % 90 degree
        oneGLRLM = rle_0(si',nl);
    case 4
        % 135 degree
        seq = zigzag(fliplr(si));
        oneGLRLM = rle_45(seq,nl);
    otherwise
        error('Only 4 directions supported')
end


% --------------------------------------------------------------------
function [I, offset, nl, gl] = ParseInputs(varargin)
% parsing parameter checking
% Inputs must be max seven item
iptchecknargin(1,7,nargin,mfilename);
%
% Check I
I = varargin{1};
iptcheckinput(I,{'logical','numeric'},{'2d','real','nonsparse'}, ...
    mfilename,'I',1);
% ------------------------
% Assign Defaults
% -------------------------
% four directions 0, 45, 90,135
offset = [1;2;3;4];
%
if islogical(I)
    nl = 2;
else
    nl = 8;
end
gl = getrangefromclass(I);

% Parse Input Arguments
if nargin ~= 1

    paramStrings = {'Offset','NumLevels','GrayLimits'};

    for k = 2:2:nargin

        param = lower(varargin{k});
        inputStr = iptcheckstrs(param, paramStrings, mfilename, 'PARAM', k);
        idx = k + 1;  %Advance index to the VALUE portion of the input.
        if idx > nargin
            eid = sprintf('Images:%s:missingParameterValue', mfilename);
            msg = sprintf('Parameter ''%s'' must be followed by a value.', inputStr);
            error(eid,'%s', msg);
        end

        switch (inputStr)

            case 'Offset'

                offset = varargin{idx};
                iptcheckinput(offset,{'logical','numeric'},...
                    {'d','nonempty','integer','real'},...
                    mfilename, 'OFFSET', idx);
                % must be row vector 
                if size(offset,2) ~= 1
                    eid = sprintf('Images:%s:invalidOffsetSize',mfilename);
                    msg = 'OFFSET must be an n x 1 array.';
                    error(eid,'%s',msg);
                end
                offset = double(offset);

            case 'NumLevels'

                nl = varargin{idx};
                iptcheckinput(nl,{'logical','numeric'},...
                    {'real','integer','nonnegative','nonempty','nonsparse'},...
                    mfilename, 'NL', idx);
                if numel(nl) > 1
                    eid = sprintf('Images:%s:invalidNumLevels',mfilename);
                    msg = 'NL cannot contain more than one element.';
                    error(eid,'%s',msg);
                elseif islogical(I) && nl ~= 2
                    eid = sprintf('Images:%s:invalidNumLevelsForBinary',mfilename);
                    msg = 'NL must be two for a binary image.';
                    error(eid,'%s',msg);
                end
                nl = double(nl);

            case 'GrayLimits'

                gl = varargin{idx};
                iptcheckinput(gl,{'logical','numeric'},{'vector','real'},...
                    mfilename, 'GL', idx);
                if isempty(gl)
                    gl = [min(I(:)) max(I(:))];
                elseif numel(gl) ~= 2
                    eid = sprintf('Images:%s:invalidGrayLimitsSize',mfilename);
                    msg = 'GL must be a two-element vector.';
                    error(eid,'%s',msg);
                end
                gl = double(gl);
        end
    end
end

function stats = grayrlprops(varargin)

%GRAYCOPROPS Properties of gray-level run-length matrix.
%  -------------------------------------------
%  STATS = GRAYCOPROPS(GLRLM,PROPERTIES) Each element in  GLRLM, (r,c),
%   is the probability occurrence of pixel having gray level values r, run-length c in the image.
%   GRAYCOPROPS is to calculate PROPERTIES.
%  -------------------------------------------
%  Requirements:
%  -------------------------------------------
%   GLRLM mustbe an cell array of valid gray-level run-length
%   matrices.Recall that a valid glrlm must be logical or numerical.
%  -------------------------------------------
%  Current supported statistics include:
%  -------------------------------------------
%   Short Run Emphasis (SRE)
%   Long Run Emphasis (LRE)
%   Gray-Level Nonuniformity (GLN)
%   Run Length Nonuniformity (RLN)
%   Run Percentage (RP)
%   Low Gray-Level Run Emphasis (LGRE)
%   High Gray-Level Run Emphasis (HGRE)
%   Short Run Low Gray-Level Emphasis (SRLGE)
%   Short Run High Gray-Level Emphasis (SRHGE)
%   Long Run Low Gray-Level Emphasis (LRLGE)
%   Long Run High Gray-Level Emphasis (LRHGE)
%  --------------------------------------------
%  Reference:
%  --------------------------------------------
%   Xiaoou Tang,Texture Information in Run-Length Matrices
%   IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL.7, NO.11,NOVEMBER 1998
% ---------------------------------------------
%  See also GRAYRLMATRIX.
% ---------------------------------------------
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
% ---------------------------------------------
% History:
% ---------------------------------------------
% Creation: beta         Date: 01/10/2007
% Revision: 1.0          Date: 12/11/2007
% 1.Accept cell input now
% 2.Using MATLAB file style
% 3.Fully vectorized programming
% 4.Fully support the IEEE reference
% 5. ...


% Check GLRLM
[GLRLM numGLRLM] = ParseInputs2(varargin{:});

% Initialize output stats structure.
% 11 statistics for each GLRLM
numStats = 11;

% % count number of GLRLM
% numGLRLM = length(GLRLM);

% Initialization default 4*11 matrix
% % stats = zeros(numGLRLM,numStats);
stats = table;

for p = 1 : numGLRLM
    %N-D indexing not allowed for sparse.

    if numGLRLM ~= 1
        % transfer to double matrix
        tGLRLM = GLRLM{p};
    else
        tGLRLM = GLRLM;
    end
    %     if numGLRLM ~= 1
    %         % transfer to double matrix
    %         tGLRLM = normalizeGLRL(GLRLM{p});
    %     else
    %         tGLRLM = normalizeGLRL(GLRLM);
    %     end
    % Get row and column subscripts of GLRLM.  These subscripts correspond to the
    % pixel values in the GLRLM.
    s = size(tGLRLM);
    % colum indicator
    c_vector =1:s(1);
    % row indicator
    r_vector =1:s(2);
    % matrix element indicator
    % Matrix form col and row: using meshgrid, you should transpose before using
    % i.e. if tGLRLM is m*n, then this function return c_matrix n*m,
    % r_matrix n*m.
    [c_matrix,r_matrix] = meshgrid(c_vector,r_vector);

    % Total number of runs
    N_runs = sum(sum(tGLRLM));

    % total number of elements
    N_tGLRLM = s(1)*s(2);

    %--------------------Prepare four matrix for speedup--------------
    % 1.Gray Level Run-Length Pixel Number Matrix
    %     p_p = calculate_p_p(tGLRLM,c_matrix');

    % 2.Gray-Level Run-Number Vector
    %   This vector represents the sum distribution of the number of runs
    %   with gray level i.
    p_g = sum(tGLRLM);

    % 3.Run-Length Run-Number Vector
    %   This vector represents the sum distribution of the number of runs
    %   with run length j.
    p_r = sum(tGLRLM,2)';

    % 4.Gray-Level Run-Length-One Vector
    %
    % p_o = tGLRLM(:,1); % Not used yet
    % ----------------------End four matrix---------------------------
    %
    %------------------------Statistics-------------------------------
    % 1. Short Run Emphasis (SRE)
    stats.SRE(p,1) = sum(p_r./(c_vector.^2))/N_runs;
    % 2. Long Run Emphasis (LRE)
    stats.LRE(p,1) = sum(p_r.*(c_vector.^2))/N_runs;
    % 3. Gray-Level Nonuniformity (GLN)
    stats.GLN(p,1) = sum(p_g.^2)/N_runs;
    % 4. Run Length Nonuniformity (RLN)
    stats.RLN(p,1) = sum(p_r.^2)/N_runs;
    % 5. Run Percentage (RP)
    stats.RP(p,1) = N_runs/N_tGLRLM;
    % 6. Low Gray-Level Run Emphasis (LGRE)
    stats.LGRE(p,1) = sum(p_g./(r_vector.^2))/N_runs;
    % 7. High Gray-Level Run Emphasis (HGRE)
    stats.HGRE(p,1) = sum(p_g.*r_vector.^2)/N_runs;
    % 8. Short Run Low Gray-Level Emphasis (SRLGE)
    stats.SGLGE(p,1) =calculate_SGLGE(tGLRLM,r_matrix',c_matrix',N_runs);
    % 9. Short Run High Gray-Level Emphasis (SRHGE)
    stats.SRHGE(p,1) =calculate_SRHGE(tGLRLM,r_matrix',c_matrix',N_runs);
    % 10. Long Run Low Gray-Level Emphasis (LRLGE)
    stats.LRLGE(p,1) =calculate_LRLGE(tGLRLM,r_matrix',c_matrix',N_runs);
    % 11.Long Run High Gray-Level Emphasis (LRHGE
    stats.LRHGE(p,1) =calculate_LRHGE(tGLRLM,r_matrix',c_matrix',N_runs);
    %----------------insert statistics----------------------------
% %     stats(p,:)=[SRE LRE GLN RLN  RP LGRE HGRE SGLGE SRHGE LRLGE  LRHGE ];
end % end all run length matrixs

%   ----------------------Utility functions--------------------
%-----------------------------------------------------------------------------
% function glrl = normalizeGLRL(glrl)
%
% % Normalize glcm so that sum(glcm(:)) is one.
% if any(glrl(:))
%   glrl = glrl ./ sum(glrl(:));
% end
% function p_p = calculate_p_p(GLRLM,c) % Note: currently not used
%
% % p_p(i; j) = GLRLM(i,j)*j
% % Each element of the matrix represents the number of pixels of run length
% % j and gray-level i. Compared to the original matrix, the new matrix gives
% % equal emphasis to all length of runs in an image.
%
% term1 =  c; % j index in matrix size
% term2 = GLRLM;
% p_p = term1 .* term2;
%---------------------------------
function SGLGE =calculate_SGLGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Short Run Low Gray-Level Emphasis (SRLGE):

term = tGLRLM./((r_matrix.*c_matrix).^2);
SGLGE= sum(sum(term))./N_runs;

%------------------------------------
function  SRHGE =calculate_SRHGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Short Run High Gray-Level Emphasis (SRHGE):
%
term  = tGLRLM.*(r_matrix.^2)./(c_matrix.^2);
SRHGE = sum(sum(term))/N_runs;
%------------------------------------
function   LRLGE =calculate_LRLGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Long Run Low Gray-Level Emphasis (LRLGE):
%
term  = tGLRLM.*(c_matrix.^2)./(r_matrix.^2);
LRLGE = sum(sum(term))/N_runs;
%---------------------------------------
function  LRHGE =calculate_LRHGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Long Run High Gray-Level Emphasis (LRHGE):
%
term  = tGLRLM.*(c_matrix.^2).*(r_matrix.^2);
LRHGE = sum(sum(term))/N_runs;
%----------------------------------------

%-----------------------------------------------------------------------------
function [glrlm num_glrlm] = ParseInputs2(varargin)
% check stability of inputs
%
% first receive all inputs
glrlm = varargin{:};
% get numbers total
num_glrlm=length(glrlm);
% then for each element, check its stability
for i=1:num_glrlm
    % The 'nonnan' and 'finite' attributes are not added to iptcheckinput because the
    % 'integer' attribute takes care of these requirements.
    % iptcheckinput(glrlm,{'cell'},{'real','nonnegative','integer'}, ...
    % mfilename,'GLRLM',1);
    iptcheckinput(glrlm{i},{'logical','numeric'},{'real','nonnegative','integer'},...
        mfilename,'GLRLM',1);
    % Cast GLRLM to double to avoid truncation by data type. Note that GLRLM is not an
    % image.
    if ~isa(glrlm,'double')
        glrlm{i}= double(glrlm{i});
    end
end

function oneglrlm = rle_0(si,NL)
% RLE   image gray level Run Length matrix for 0degree
%    
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
% History:
%  -------
% Creation: beta  Date: 01/11/2007 
% Revision: 1.0   Date: 10/11/2007


% Assure row number is exactly the gray level
[m,n]=size(si);

oneglrlm=zeros(NL,n);

for i=1:m
    x=si(i,:);
    % run length Encode of each vector
    index = [ find(x(1:end-1) ~= x(2:end)), length(x) ];
    len = diff([ 0 index ]); % run lengths
    val = x(index);          % run values
    temp =accumarray([val;len]',1,[NL n]);% compute current numbers (or contribution) for each bin in GLRLM
    oneglrlm = temp + oneglrlm; % accumulate each contribution
end

