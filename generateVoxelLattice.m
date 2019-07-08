function [voxel,Density] = generateVoxelLattice(n,radius,node,strut)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n is the number of voxel along each axis
% address is the file location of wireframe
% density is the relative density

% Copyright (c) 2018, Guoying Dong
% (https://github.com/GuoyingDong/homogenization)
% Modified: Yu-Chin Chan (ychan@u.northwestern.edu), 3/26/2019, 6/18/2019
%     - 3/26/2019: remove need for text file as input
%     - 6/18/2019: vectorize for significant speedup (e.g. 0.8sec vs. 9sec for n=40)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vox_size = 1/n;               % initial size of voxels
voxel = zeros(n,n,n);      % initial grid with zeros
%% generate a list of centers of voxel
voxel_c = zeros(n^3,6);   
p = 0;                    % p count the number of all voxels
for i = 1:n               % i for z axis
    for j = 1:n           % j for y axis
        for k = 1:n       % k for x axis
            p = p + 1;
            voxel_c(p,1:3) = [k,j,i];  % save index along x,y,z axis
            % save coordinate along x,y,z axis
            voxel_c(p,4:6) = [(k-0.5)*vox_size,(j-0.5)*vox_size,(i-0.5)*vox_size];         
        end
    end
end
voxel_i = sub2ind(size(voxel), voxel_c(:,1), voxel_c(:,2), voxel_c(:,3));
start_n = node(strut(:,1),:);
end_n = node(strut(:,2),:);

%% Get the voxel close the the strut witnin a certain distance
for i = 1:length(strut)
    alpha = acosd( sum((voxel_c(:,4:6) - start_n(i,:)) .* (end_n(i,:) - start_n(i,:)), 2)...
        ./ (vecNorm(voxel_c(:,4:6) - start_n(i,:)) .* vecNorm(end_n(i,:) - start_n(i,:))) );
    beta  = acosd( sum((voxel_c(:,4:6) - end_n(i,:)) .* (start_n(i,:) - end_n(i,:)), 2)...
        ./ (vecNorm(voxel_c(:,4:6) - end_n(i,:)) .* vecNorm(start_n(i,:) - end_n(i,:))) );

    % if it is acute angle, distance to node
    distance = min(vecNorm(voxel_c(:,4:6) - start_n(i,:)),...
        vecNorm(voxel_c(:,4:6) - end_n(i,:)));
    % if not acute angle, distance to line
    obtuse = (alpha<90 & beta<90);
    temp = vecNorm( ...
        cross(repmat(end_n(i,:) - start_n(i,:), p, 1), voxel_c(:,4:6) - start_n(i,:), 2) )...
        ./ vecNorm(end_n(i,:) - start_n(i,:));
    distance(obtuse) = temp(obtuse);

    % if distance less than radius, activate it
    temp = zeros(p,1);
    active = (distance <= radius);
    temp(active) = 1;
    temp_voxel = zeros(size(voxel));
    temp_voxel(voxel_i) = temp;
    voxel = temp_voxel | voxel;
end
Density = sum(sum(sum(voxel)))/numel(voxel); % calculate the relative density
end

function new_norm = vecNorm(A)
new_norm = sqrt(sum(A.^2, 2));
end