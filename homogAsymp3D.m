function CH = homogAsymp3D(ll, vox, props0, def, solver)
% =========================================================================
% Asymptotic homogenization of 2-phase 3D static linear elastic periodic 
% structures using eight-node hexahedral elements
% 
% Note: only for 2-phases with same Poisson's ratio when using def='youngs'
%       not an issue with def='lame' option
% 
% Syntax: CH = homogAsymp3D([1,1,1], vox, [1e-9,params.E; params.v,params.v], 'youngs', 'pcg');
% 
% Inputs:
%     ll     - lengths of sides of unit cell, [lx,ly,lz]
%     voxel  - binary or 0/1 3D matrix (0 = material 1, 1 = material 2)
%     def    - definition of isotropic constituent material properties:
%               1) 'lame' to use Lame's parameters lambda and mu
%               2) 'youngs' to use Young's modulus E and Poisson's ratio nu
%     props0 - if opt = 'lame', props0 = [lambda1, lambda2; mu1, mu2];
%              if opt = 'youngs', props0 = [E1, E2; nu1, nu2];
%     solver - matrix solver option:
%               1) 'pcg' to use preconditioned conjugate gradient solver
%                  (recommended for large number of elements)
%               2) 'direct' to use direct matrix inverse solver
% 
% Author: Yu-Chin Chan (ychan@u.northwestern.edu), 4/12/2019
% Last updated: 4/25/2019
% 
% Based on code by:
% [1] Dong G, Tang Y, Zhao Y, "A 149 Line Homogenization Code for 
%     Three-Dimensional Cellular Materials Written in matlab",
%     J Eng Mater Technol, 2018, doi:10.1115/1.4040555,
%     https://github.com/GuoyingDong/homogenization.
% =========================================================================
%% INITIALIZE
if ~exist('def', 'var')
    def = 'youngs';
end
if ~exist('solver', 'var')
    solver = 'pcg';
end
[nelx, nely, nelz] = size(vox); %size of voxel model along x,y and z axis
dx = ll(1)/nelx; dy = ll(2)/nely; dz = ll(3)/nelz;
nel = nelx*nely*nelz;
% Node numbers and element degrees of freedom for full (not periodic) mesh
nodenrs = reshape(1:(1+nelx)*(1+nely)*(1+nelz),1+nelx,1+nely,1+nelz);
edofVec = reshape(3*nodenrs(1:end-1,1:end-1,1:end-1)+1,nel,1);
addx = [0 1 2 3*nelx + [3 4 5 0 1 2] -3 -2 -1];
addxy = 3*(nely+1)*(nelx+1) + addx;
edofMat = repmat(edofVec,1,24) + repmat([addx addxy],nel,1);

%% IMPOSE PERIODIC BOUNDARY CONDITIONS
% Use original edofMat to index into list with the periodic dofs
nn = (nelx+1)*(nely+1)*(nelz+1); % Total number of nodes
nnP = (nelx)*(nely)*(nelz);      % Total number of unique nodes
nnPArray = reshape(1:nnP, nelx, nely, nelz);
% Extend with a mirror of the back border
nnPArray(end+1,:,:) = nnPArray(1,:,:);
% Extend with a mirror of the left border
nnPArray(:, end+1, :) = nnPArray(:,1,:);
% Extend with a mirror of the top border
nnPArray(:, :, end+1) = nnPArray(:,:,1);
% Make a vector into which we can index using edofMat:
dofVector = zeros(3*nn, 1);
dofVector(1:3:end) = 3*nnPArray(:)-2;
dofVector(2:3:end) = 3*nnPArray(:)-1;
dofVector(3:3:end) = 3*nnPArray(:);
edof = dofVector(edofMat);
ndof = 3*nnP;

%% ASSEMBLE GLOBAL STIFFNESS MATRIX K AND LOAD VECTOR F
% Indexing vectors
iK = kron(edof,ones(24,1))';
jK = kron(edof,ones(1,24))';
iF = repmat(edof',6,1);
jF = [ones(24,nel); 2*ones(24,nel); 3*ones(24,nel);...
    4*ones(24,nel); 5*ones(24,nel); 6*ones(24,nel);];
% Assemble stiffness matrix and load vector
if strcmp(def, 'lame')
    % Material properties assigned to voxels with materials
    lambda = props0(1,:); mu = props0(2,:);
    lambda = lambda(1)*(vox==0) + lambda(2)*(vox==1);
    mu = mu(1)*(vox==0) + mu(2)*(vox==1);
    % Unit element stiffness matrix and load
    [keLambda, keMu, feLambda, feMu] = assemble_lame(dx/2, dy/2, dz/2);
    ke = keMu + keLambda; % Here the exact ratio does not matter, because
    fe = feMu + feLambda; % it is reflected in the load vector
    sK = keLambda(:)*lambda(:).' + keMu(:)*mu(:).';
    sF = feLambda(:)*lambda(:).' + feMu(:)*mu(:).';
elseif strcmp(def, 'youngs')
    E = props0(1,:); E = E(1)+vox.*(E(2)-E(1)); % SIMP
    nu = props0(2,:);
    % Unit element stiffness matrix and load
    [ke, fe] = assemble_youngs(nu, dx/2, dy/2, dz/2);
    sK = ke(:)*E(:)';
    sF = fe(:)*E(:)';
else
    error('unavailable option for constituent properties definition')
end
% Global stiffness matrix
K = sparse(iK(:), jK(:), sK(:), ndof, ndof);
K = (K+K')/2;
% Six load cases corresponding to the six strain cases
F  = sparse(iF(:), jF(:), sF(:), ndof, 6);

%% SOLUTION
activedofs = edof((vox==0 | vox==1),:);
activedofs = sort(unique(activedofs(:)));
X = zeros(ndof,6);
if strcmp(solver, 'pcg')
    % solve using PCG method, remember to constrain one node
    L = ichol(K(activedofs(4:end),activedofs(4:end))); % preconditioner
    for i = 1:6 % run once for each loading condition
        [X(activedofs(4:end),i),~,~,~] = pcg(K(activedofs(4:end),...
        activedofs(4:end)),F(activedofs(4:end),i),1e-10,300,L,L');
    end
elseif strcmp(solver, 'direct')
    % solve using direct method
    X(activedofs(4:end),:) = K(activedofs(4:end),activedofs(4:end))...
        \F(activedofs(4:end),:);
else
    error('unavailable option for solver')
end

%% ASYMPTOTIC HOMOGENIZATION
% The displacement vectors corresponding to the unit strain cases
X0 = zeros(nel, 24, 6);
% The element displacements for the six unit strains
X0_e = zeros(24, 6);
% fix degrees of nodes [1 2 3 5 6 12];
X0_e([4 7:11 13:24],:) = ke([4 7:11 13:24],[4 7:11 13:24])...
                           \fe([4 7:11 13:24],:);
X0(:,:,1) = kron(X0_e(:,1)', ones(nel,1)); % epsilon0_11 = (1,0,0,0,0,0)
X0(:,:,2) = kron(X0_e(:,2)', ones(nel,1)); % epsilon0_22 = (0,1,0,0,0,0)
X0(:,:,3) = kron(X0_e(:,3)', ones(nel,1)); % epsilon0_33 = (0,0,1,0,0,0)
X0(:,:,4) = kron(X0_e(:,4)', ones(nel,1)); % epsilon0_12 = (0,0,0,1,0,0)
X0(:,:,5) = kron(X0_e(:,5)', ones(nel,1)); % epsilon0_23 = (0,0,0,0,1,0)
X0(:,:,6) = kron(X0_e(:,6)', ones(nel,1)); % epsilon0_13 = (0,0,0,0,0,1)
CH = zeros(6);
volume = prod(ll);
% Homogenized elasticity tensor
if strcmp(def, 'lame')
    for i = 1:6
        for j = 1:6
            sum_L = ((X0(:,:,i) - X(edof+(i-1)*ndof))*keLambda).*...
                (X0(:,:,j) - X(edof+(j-1)*ndof));
            sum_M = ((X0(:,:,i) - X(edof+(i-1)*ndof))*keMu).*...
                (X0(:,:,j) - X(edof+(j-1)*ndof));
            sum_L = reshape(sum(sum_L,2), nelx, nely, nelz);
            sum_M = reshape(sum(sum_M,2), nelx, nely, nelz);
            CH(i,j) = 1/volume*sum(sum(sum(lambda.*sum_L + mu.*sum_M)));
        end
    end
elseif strcmp(def, 'youngs')
    for i = 1:6
        for j = 1:6
            sum_XkX = ((X0(:,:,i) - X(edof+(i-1)*ndof))*ke).*...
                (X0(:,:,j) - X(edof+(j-1)*ndof));
            sum_XkX = reshape(sum(sum_XkX,2), nelx, nely, nelz);
            CH(i,j) = 1/volume*sum(sum(sum(sum_XkX.*E)));
        end
    end
end
end

%% COMPUTE UNIT ELEMENT STIFFNESS MATRIX AND LOAD VECTOR
function [keLambda, keMu, feLambda, feMu] = assemble_lame(a, b, c)
% Initialize
keLambda = zeros(24,24); keMu = zeros(24,24);
feLambda = zeros(24,6); feMu = zeros(24,6);
ww = [5/9, 8/9, 5/9];
J_ = [-a a a -a -a a a -a; -b -b b b -b -b b b; -c -c -c -c c c c c]';
% Constitutive matrix contributions
CMu = diag([2 2 2 1 1 1]); CLambda = zeros(6); CLambda(1:3,1:3) = 1;
% Three Gauss points in both directions
xx = [-sqrt(3/5), 0, sqrt(3/5)]; yy = xx; zz = xx;
for ii = 1:length(xx)
    for jj = 1:length(yy)
        for kk = 1:length(zz)
            % integration point
            x = xx(ii); y = yy(jj); z = zz(kk);
            % stress strain displacement matrix
            [B, J] = strain_disp_matrix(x, y, z, J_);
            % Weight factor at this point
            weight = det(J) * ww(ii) * ww(jj) * ww(kk);
            % Element matrices
            keLambda = keLambda + weight * B' * CLambda * B;
            keMu = keMu + weight * B' * CMu * B;
            % Element loads
            feLambda = feLambda + weight * B' * CLambda;       
            feMu = feMu + weight * B' * CMu; 
        end
    end
end
end

function [ke, fe] = assemble_youngs(nu, a, b, c)
%  Initialize
ww = [5/9, 8/9, 5/9];
J_ = [-a a a -a -a a a -a; -b -b b b -b -b b b; -c -c -c -c c c c c]';
ke = zeros(24,24); fe = zeros(24,6);
% Constitutive matrix with unit Young's modulus
nu = nu(2); %TODO multi-material nu
C = diag([nu nu 0 0 0],1) + diag([nu 0 0 0],2); C = C+C';
C = C + diag([repmat(1-nu,3,1); repmat((1-2*nu)/2,3,1)]);
C = C / ((1+nu)*(1-2*nu));
% Three Gauss points in both directions
xx = [-sqrt(3/5), 0, sqrt(3/5)]; yy = xx; zz = xx;
for ii = 1:length(xx)
    for jj = 1:length(yy)
        for kk = 1:length(zz)
            % integration point
            x = xx(ii); y = yy(jj); z = zz(kk);
            % stress strain displacement matrix
            [B, J] = strain_disp_matrix(x, y, z, J_);
            % Weight factor at this point
            weight = det(J) * ww(ii) * ww(jj) * ww(kk);
            % Element matrices
            ke = ke + weight * B' * C * B;
            % Element loads
            fe = fe + weight * B' * C;       
        end
    end
end
end

function [B, J] = strain_disp_matrix(x, y, z, J_)
%stress strain displacement matrix
qx = [ -((y-1)*(z-1))/8, ((y-1)*(z-1))/8, -((y+1)*(z-1))/8,...
    ((y+1)*(z-1))/8, ((y-1)*(z+1))/8, -((y-1)*(z+1))/8,...
    ((y+1)*(z+1))/8, -((y+1)*(z+1))/8];
qy = [ -((x-1)*(z-1))/8, ((x+1)*(z-1))/8, -((x+1)*(z-1))/8,...
    ((x-1)*(z-1))/8, ((x-1)*(z+1))/8, -((x+1)*(z+1))/8,...
    ((x+1)*(z+1))/8, -((x-1)*(z+1))/8];
qz = [ -((x-1)*(y-1))/8, ((x+1)*(y-1))/8, -((x+1)*(y+1))/8,...
    ((x-1)*(y+1))/8, ((x-1)*(y-1))/8, -((x+1)*(y-1))/8,...
    ((x+1)*(y+1))/8, -((x-1)*(y+1))/8];
J = [qx; qy; qz] * J_; % Jacobian
qxyz = J \ [qx; qy; qz];
B_e = zeros(6,3,8);
for i_B = 1:8
    B_e(:,:,i_B) = [qxyz(1,i_B)   0             0;
                    0             qxyz(2,i_B)   0;
                    0             0             qxyz(3,i_B);
                    qxyz(2,i_B)   qxyz(1,i_B)   0;
                    0             qxyz(3,i_B)   qxyz(2,i_B);
                    qxyz(3,i_B)   0             qxyz(1,i_B)];
end
B = [B_e(:,:,1) B_e(:,:,2) B_e(:,:,3) B_e(:,:,4) B_e(:,:,5)...
    B_e(:,:,6) B_e(:,:,7) B_e(:,:,8)];
end