% Yu-Chin Chan, (ychan@u.northwestern.edu), 7/8/2019

clearvars, close all, clc

% define GRID (nx3 array of node coordinates) and STRUT (sx2 array of node connectivity)
load('topology/grid_octet_skel.mat');
% voxelize
res = 40; % number of voxels per side
rad = 0.1; % radius of struts
tic
[vox, dens] = generateVoxelLattice(res, rad, GRID, STRUT);
toc

% alternatively, define voxels directly
% load('topology/grid_octet_vox.mat');
% dens = sum(sum(sum(vox)))/numel(vox);

% lengths of sides of unit cell
ll = [1,1,1];

% properties of isotropic constituent material properties
E = [1e-9, 2e9]; % E1, E2
nu = [0.33, 0.33]; % nu1, nu2
lam = nu.*E ./ ((1+nu).*(1-2*nu));
mu = E ./ (2*(1+nu));

% two options to define constituent materials: 'young's or 'lame'
% changes how stiffness matrix is assembled.
def = 'youngs'; props0 = [E; nu];   % with young's modulus and poisson's ratio
% def = 'lame'; props0 = [lam; mu]; % with lame's parameters

% two options for solver: 'pcg' or 'direct'
solver = 'pcg';

% two options to format effective property results: 'vec' or 'struct'
outOption = 'struct';

% options to print results and/or plot Young's modulus surface
dispFlag = 1;
plotFlag = 1;

% run homogenization
tic
CH = homogAsymp3D(ll, vox, props0, def, solver);
[props, SH] = evaluateCH(CH, dens, outOption, dispFlag);
toc
if plotFlag
    visual(CH);
    axis equal; colorbar;
end