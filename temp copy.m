% all stats tests
%% init vars + names
% imgs = {'baselineneigbdfltext'};
% imgs = {'texturedflair'};

imgs = {'baselinet', 'neighbaseline', 'texturedflair', 'baselineneigbdfltext'};
% 15 files for each, generate 5 values for each file. 5*15*4 = 300
patients = [1,2,4,5,6,7,8,10,12,14,15,16,17,18,19]; %no AUC = 0 no disease progression from combine.m
% patients = [1];
AUC = [];
PR = [];
PPV = [];
DSI = [];
SE = [];
OPT = [];
%% calculate values for each
for i = 1:length(imgs)
    auc = [];
    pr = [];
    ppv = [];
    dsi = [];
    se = [];
    opt = [];
    name = imgs{i};
    for J = patients %J is patient number
        %% load actual predictions
        num = num2str(J);
        folder = strcat('../genMSdata2/patient', num);
        filePattern = fullfile(folder, strcat('*', name,'*')); % Change to whatever pattern you need.
        theFiles = dir(filePattern);
        s.patient_number = num;
        fullFileName = [];
        for k = 1 : length(theFiles)
            baseFile = theFiles(k).name;
            if k == 1
                a = getfield(load_nii(fullfile(folder, baseFile)),'img');
            elseif k == 2
                b = getfield(load_nii(fullfile(folder, baseFile)),'img');
            elseif k == 3
                c = getfield(load_nii(fullfile(folder, baseFile)),'img');
            end
        end
        gt3 = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_gt3.nii')),'img')==1;
        if strfind(name, 'text')
            mask = zeros(size(gt3));
            directoryName = strcat('texturepatches3/patient', num);
            inds = h5read(strcat(directoryName, '/patient', num,'_unbalancedtextureProps.h5'),'/inds');
            mask(inds)= 1;
            mask = logical(mask);
      else
%             fprintf('yay')
            mask = getfield(load_nii(strcat('../MSpatientdata/patient',num,'/patient',num,'_brainmask.nii.gz')),'img')>0;
      
        end
        gt = gt3(mask);
        a = a(mask);
        b = b(mask);
        c = c(mask);
    
    iter = {a, b, c};
        for i = 1:length(iter) 
           new = iter{i};
           new(isnan(new) |isinf(new))=0;
             %% AUC, OPT
            [X,Y,T,AUCROC,OPTROCPT] = perfcurve(gt, new, true,'xCrit','FPR','TVals',min(new):range(new)/100:max(new));%testGroundTruth{J}
            AUCROC;
            auc = [auc; AUCROC];
            thr = T(find(X>=OPTROCPT(1),1));
            opt = [opt; thr];
            plot(X,Y);
            axis square
%%         AUCPR         
%             [X,Y,T,AUCPR,OPTROCPT] = perfcurve(gt, new, true,'xCrit','PPV');%testGroundTruth{J
%             pr = [pr; AUCPR];
%                 hold on;
%                 plot(X,Y);
%             %     title = strcat
%                 title('Test Data PR');
%                 axis square
%             %% dsi
%                 dice = @(auto_seg,manual_seg) 2*nnz(auto_seg & manual_seg)/(nnz(auto_seg) + nnz(manual_seg));
%                sens = @(auto_seg,manual_seg) nnz(auto_seg & manual_seg)/nnz(manual_seg); %here tp/fn + tp = recall
%                 ppval = @(auto_seg,manual_seg) nnz(auto_seg & manual_seg)/nnz(auto_seg); %here tp/fp + tp = precision
% 
%                 dsitest = new;
%                 dsitest(dsitest > 0.75) = 1;
%                 dsitest(dsitest < 0.75) = 0;
%                 dsi =[dsi; dice(dsitest, gt)];
%                 se =[se; sens(dsitest, gt)];
%                 ppv =[ppv; ppval(dsitest, gt)];


                
        end
    end
AUC = [AUC, auc];
PR = [PR, pr];
PPV = [PPV, ppv];
DSI = [DSI, dsi];
SE = [SE, se];
OPT = [OPT, opt];
end

% [p,h] = signrank(testAUCs_baseline,testAUCs)