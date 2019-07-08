# Asymptotic Homogenization:

**Original code by:**

Please cite the following article if you use this code in your publications:

Dong G, Tang Y, Zhao Y. A 149 Line Homogenization Code for Three-Dimensional Cellular Materials Written in matlab. ASME. J. Eng. Mater. Technol. 2018;141(1):011005-011005-11. doi:10.1115/1.4040555.

Homogenization code for 3D lattice structure. For more information you can visit http://www.intralatticepro.com/

![alt text](https://github.com/lychan110/homogAsymp3D/blob/master/image/homogenization.JPG)

**Changes to original codes:**

`homo3d.m` is now `homogAsymp3D.m`: calculate the homogenized material property of 3d lattice structures.
* Improvements:
    1. Option added to assemble the stiffness matrix using either Young's modulus/Poisson's ratio or Lame's parameters.
    2. Option added to more easily switch between pcg or direct matrix solver.
* Note: this code can be used for any voxel structure, not just lattice structures. To voxelize mesh structures in Matlab, I recommend [this code by Adam Aitkenhead](https://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation).

`GenerateVoxel.m` is now `generateVoxelLattice.m`: generate the voxel model of lattice structures for `homogAsymp3D.m`.
* Improvements:
    1. Instead of a text file, the inputs are now arrays for nodes and struts. This negates the need to read/write files during iterative optimization.
    2. The voxel activation process has been vectorized for a significant speedup. For the included grid octet example, the vectorized code takes around 0.8 sec while the original code needs around 9 sec.

`visual.m`: plot the 3D Young's modulus surface indicating E along all directions.
* No change.

# Example:

Example using `youngs` and `pcg` options (see `runHomogAsymp3D.m` for full example with other options):

The structure:

<img src="https://github.com/lychan110/homogAsymp3D/blob/master/image/grid_octet.png" width="200">

```MATLAB
% lengths of sides of unit cell
ll = [1,1,1];

% properties of isotropic constituent material properties
E = [1e-9, 2e9]; % E1, E2
nu = [0.33, 0.33]; % nu1, nu2
lam = nu.*E ./ ((1+nu).*(1-2*nu));
mu = E ./ (2*(1+nu));

% options for homogenization
def = 'youngs'; props0 = [E; nu];   % with young's modulus and poisson's ratio
solver = 'pcg';
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
    axis equal
end
```

The result:

```
--------------------------EFFECTIVE PROPERTIES--------------------------
Density: 0.4407
Youngs Modulus:____E11_____|____E22_____|____E33_____
                3.2813e+08 | 3.2813e+08 | 3.2813e+08

Shear Modulus:_____G23_____|____G31_____|____G12_____
                1.3790e+08 | 1.3790e+08 | 1.3790e+08

Poissons Ratio:____v12_____|____v13_____|____v23_____
                    0.2691 |     0.2691 |     0.2691
               ____v21_____|____v31_____|____v32_____
                    0.2691 |     0.2691 |     0.2691
------------------------------------------------------------------------    
```

The Young's modulus surface:

<img src="https://github.com/lychan110/homogAsymp3D/blob/master/image/grid_octet_surface.png" width="400">
