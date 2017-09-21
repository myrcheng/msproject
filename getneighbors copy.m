function [ret] = getneighbors(input)
size(input)    
global neighbors;
    global index;
    global neighbormask;
    index = index+1;
    temp = input;
%     size(temp)
%     maskTemp = input .* neighbormask;
%% % delete column if 5th value = 0 (get rid of voxels not in the mask)
%     [r, c]=find(maskTemp(5,:)==0); 
%     temp(:,c)=[];
    %% 
%     temp(5,:) = []; %delete duplicate center values
%     [x, y, z] = ind2sub(size(neighbormask), index);
    selected_voxels = find(neighbormask);
    if(find(selected_voxels == index))
        neighbors = [neighbors, temp];
    end
    %% 
        ret = ones(1, length(input));

end
