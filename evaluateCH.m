function [props, SH] = evaluateCH(CH, dens, outOption, dispFlag)
% =========================================================================
% Calculate mechanical moduli (elastic, shear, Poisson's ratio) 
% when given orthotropic elasticity tensor.
% 
% Inputs: CH        - elasticity tensor
%         dens      - density of structure
%         outOption - type of output:
%                         outOption = 'struct' to return properties as structure
%                         outOption = 'vec' to return as row vector
%         dispFlag  - option to print out properties to command window
%                         dispFlag = 1 to print
%                         dispFlag = 0 to not print
% 
% Output: props     - if outOption = 'struct', props is structure containing all info:
%                         CH - elasticity/stiffness tensor
%                         SH - compliance tensor
%                         EH = [E1, E2, E3]; % Young's modulus
%                         GH = [G23, G31, G12]; % shear modulus
%                         vH = [v12, v13, v23; v21, v31, v32]; % Poisson's ratio
%                   - if outOption = 'vec', props = [EH, GH, vH(:)', dens];
%         SH        - compliance tensor, optional output
%                         useful if using 'vec', which doesn't save SH to props
% 
% Author: Yu-Chin Chan (ychan@u.northwestern.edu), 4/26/2019
% Last updated: 5/15/2019, 6/18/2019
% =========================================================================
[U,S,V] = svd(CH);
sigma = diag(S);
k = sum(sigma > 1e-15);
SH = (U(:,1:k) * diag(1./sigma(1:k)) * V(:,1:k)')'; % more stable SVD (pseudo)inverse
EH = [1/SH(1,1), 1/SH(2,2), 1/SH(3,3)]; % effective Young's modulus
GH = [1/SH(4,4), 1/SH(5,5), 1/SH(6,6)]; % effective shear modulus
vH = [-SH(2,1)/SH(1,1), -SH(3,1)/SH(1,1), -SH(3,2)/SH(2,2);
      -SH(1,2)/SH(2,2), -SH(1,3)/SH(3,3), -SH(2,3)/SH(3,3)]; % effective Poisson's ratio
if strcmp(outOption, 'struct')
    props = struct('CH',CH, 'SH',SH, 'EH',EH, 'GH',GH, 'vH',vH, 'density',dens);
elseif strcmp(outOption, 'vec')
    props =  [EH, GH, vH(:)', dens];
end
if dispFlag
    fprintf('\n--------------------------EFFECTIVE PROPERTIES--------------------------\n')
    fprintf('Density: %6.4f\n', dens)
    fprintf('Youngs Modulus:____E11_____|____E22_____|____E33_____\n')
    fprintf('                %7.4e | %7.4e | %7.4e\n\n', EH)
    fprintf('Shear Modulus:_____G23_____|____G31_____|____G12_____\n')
    fprintf('                %7.4e | %7.4e | %7.4e\n\n', GH)
    fprintf('Poissons Ratio:____v12_____|____v13_____|____v23_____\n')
    fprintf('                %10.4f | %10.4f | %10.4f\n', vH(1,:))
    fprintf('               ____v21_____|____v31_____|____v32_____\n')
    fprintf('                %10.4f | %10.4f | %10.4f\n', vH(2,:))
    disp('------------------------------------------------------------------------')
end
end