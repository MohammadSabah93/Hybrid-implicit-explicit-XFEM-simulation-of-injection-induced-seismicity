% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function [Sxx,Sxy,Syy,Svm, S1, S2, PP] = elemStress(params, DISPLACEMENT, connectivity, node_coordinates, PSI, NODES, CRACK, omega, domain_length, nex, elementType)
% This function calculates the stress distribution within each element from
% the nodal displacements.  The stresses in the xx, yy and xy directions 
% are calculated.

% Extract material properties
Yu = params.Yu;  % Young's modulus
nu = params.nu;  % Poisson's ratio
Sxxr=params.Sxxr;
Syyr=params.Syyr;
Sxyr=params.Sxyr;

lXElem = domain_length/nex; 
nElements = size(connectivity,1);
nNodes = size(node_coordinates,1);
D = (Yu / ((1 + nu) * (1 - 2 * nu))) * ...
                       [1 - nu, nu, 0; 
                       nu, 1 - nu, 0; 
                       0, 0, (1 - 2 * nu) / 2];


m = size(CRACK,1);                                                          % Determine number of data points defining crack
nCT = 2;
if m > 0
    if nCT == 1
        xCT = CRACK(m,1);                                                   % X-coordinate of crack tip
        yCT = CRACK(m,2);                                                   % Y-coordinate of crack tip
    elseif nCT == 2
        xCT = [CRACK(m,1) CRACK(1,1)];                                      % X-coordinates of crack tips
        yCT = [CRACK(m,2) CRACK(1,2)];                                      % Y-coordinates of crack tips
    end
end

stressXX = zeros(nElements,4);                                                  % Create a matrix for storing nodal stress values
stressXY = zeros(nElements,4);                                                  % Create a matrix for storing nodal stress values
stressYY = zeros(nElements,4);                                                  % Create a matrix for storing nodal stress values
stressVM = zeros(nElements,4);                                                  % Create a matrix for storing nodal stress values

gp = [-1 -1; 1 -1; 1  1;-1  1];                                             % Gauss points defining the nodes of the elements

% Calculate stress at each gauss point
for iElem = 1:nElements

elementNodes = connectivity(iElem, :)';
nodesPerElement =size(connectivity(iElem, :),2);
Xe = [node_coordinates(connectivity(iElem, :), 1), node_coordinates(connectivity(iElem, :), 2)];
CTN = nnz(NODES(elementNodes,4));                                   % Number of nodes with crack tip enrichment    
HEN = nnz(NODES(elementNodes,2));                                   % Number of nodes with Heaviside enrichment
NEN = HEN+CTN;                                                      % Number of enriched nodes

local  = reshape([2 * connectivity(iElem, :) - 1; 2 * connectivity(iElem, :)], [], 1);

% Elemental displacement = [u1;v1;u2;v2;u3;v3;u4;v4]
U = DISPLACEMENT(local);

    if (NEN == 0)                                                           % Unenriched element
        for iGP = 1:length(gp)

            xi = gp(iGP,1); eta = gp(iGP,2);                                % Current gauss point
            [~, dNLocal]=ShapeFunction([xi,eta],elementType);
            jacobian = dNLocal * Xe;
            dNGlobal = jacobian \ dNLocal;
            
            Bu = zeros(3, 2 * nodesPerElement);
            Bu(1, 1:2:end) = dNGlobal(1, :);  % ∂N/∂x for x-direction
            Bu(2, 2:2:end) = dNGlobal(2, :);  % ∂N/∂y for y-direction
            Bu(3, 1:2:end) = dNGlobal(2, :);  % ∂N/∂y for x-direction (shear)
            Bu(3, 2:2:end) = dNGlobal(1, :);  % ∂N/∂x for y-direction (shear)
            
            stress = D*Bu*U;
            stressXX(iElem,iGP) = stress(1);
            stressYY(iElem,iGP) = stress(2);
            stressXY(iElem,iGP) = stress(3);
            stressVM(iElem,iGP) = sqrt(stress(1)^2+stress(2)^2-stress(1)*stress(2)+3*stress(3)^2);            
        end
    else
        H_nodes = heavisideFunction(PSI(elementNodes), 'type2'); 
        Uenr = [];
        for iGP = 1:length(gp)
            xi = gp(iGP,1); eta = gp(iGP,2);

           [N, dNLocal]=ShapeFunction([xi,eta],elementType);

            Xgp = N * Xe(:,1);                          % The global X for the current gauss point
            Ygp = N * Xe(:,2);                          % The global Y for the current gauss point
            jacobian = dNLocal * Xe;
            dNGlobal = jacobian \ dNLocal;
            
            Bu = zeros(3, 2 * nodesPerElement);
            Bu(1, 1:2:end) = dNGlobal(1, :);  % ∂N/∂x for x-direction
            Bu(2, 2:2:end) = dNGlobal(2, :);  % ∂N/∂y for y-direction
            Bu(3, 1:2:end) = dNGlobal(2, :);  % ∂N/∂y for x-direction (shear)
            Bu(3, 2:2:end) = dNGlobal(1, :);  % ∂N/∂x for y-direction (shear)

            iB = 1; iLoc = 1; Benr = []; 
            for iN = 1:nodesPerElement
                 if NODES(elementNodes(iN),2) ~= 0
                    psi1 = PSI(elementNodes(1));                                         % Psi level set value at node 1
                    psi2 = PSI(elementNodes(2));                                         % Psi level set value at node 2
                    psi3 = PSI(elementNodes(3));                                         % Psi level set value at node 3
                    psi4 = PSI(elementNodes(4));                                         % Psi level set value at node 4
                    psi  = N(1)*psi1+N(2)*psi2+N(3)*psi3+N(4)*psi4;         % Psi level set value at current gauss point
                    if psi == 1e-6, psi = 0; end

                    Hgp  = sign(psi);
                    if Hgp == 0, Hgp = sign(nonzeros([psi1 psi2 psi3 psi4])); end   
                    if length(Hgp) > 1, Hgp = sign(Hgp(1)); end 
                    
                    H   = (Hgp-H_nodes(iN));                                       % Shifted Heaviside value

                    Ba = [dNGlobal(1,iN)*H  0;
                          0    dNGlobal(2,iN)*H ;
                          dNGlobal(2,iN)*H dNGlobal(1,iN)*H];
                    
                    Benr(:,iB:(iB+1)) = Ba;
                    iB = iB+2;

                    if iGP == 1
                        Uenr(iLoc:(iLoc+1)) = [DISPLACEMENT(2*NODES(elementNodes(iN),2)-1,1) DISPLACEMENT(2*NODES(elementNodes(iN),2),1)];
                        iLoc = iLoc+2;
                    end
                 end
            end
                 for iN = 1:nodesPerElement
                   if NODES(elementNodes(iN),4) ~= 0
                   if nCT == 1
                        X     = Xgp-xCT;                                    % Horizontal distance from crack tip to gauss point
                        Y     = Ygp-yCT;                                    % Vertical distance from crack tip to gauss point
                        CCS   = [cos(omega) sin(omega);-sin(omega) cos(omega)];
                        XYloc = CCS*[X Y]';                                 % Change to crack tip coordinates
                        r     = sqrt(XYloc(1)^2+XYloc(2)^2);                % Radius from crack tip to current gauss point
                        if r < 0.001*lXElem; r = 0.05*lXElem; end
                        theta = atan2(XYloc(2),XYloc(1));                   % Angle from crack tip to current gauss point
                    elseif nCT == 2
                        X1  = Xgp-xCT(1);
                        Y1  = Ygp-yCT(1);
                        X2  = Xgp-xCT(2);
                        Y2  = Ygp-yCT(2);
                        CCS = [cos(omega(1)) sin(omega(1));-sin(omega(1)) cos(omega(1))];
                        XY1 = CCS*[X1 Y1]';
                        CCS = [cos(omega(2)) sin(omega(2));-sin(omega(2)) cos(omega(2))];
                        XY2 = CCS*[X2 Y2]';
                        r1  = sqrt(XY1(1)^2+XY1(2)^2);                      % Radius from crack tip to current gauss point
                        r2  = sqrt(XY2(1)^2+XY2(2)^2);
                        if r1 > r2
                            r = r2; theta = atan2(XY2(2),XY2(1));
                            CCS = [cos(omega(2)) sin(omega(2));-sin(omega(2)) cos(omega(2))];                            
                        elseif r2 > r1
                            r = r1; theta = atan2(XY1(2),XY1(1));
                            CCS = [cos(omega(1)) sin(omega(1));-sin(omega(1)) cos(omega(1))];
                        end
                        if r < 0.001*lXElem; r = 0.05*lXElem; end     
                    end

              c = 1/(2*sqrt(r)); ct = CCS(1,1); st = CCS(1,2);         % Constants
            

            a1gp = sqrt(r)*sin(theta/2);                        % Node 1 crack tip enrichment value
            a2gp = sqrt(r)*cos(theta/2);                        % Node 2 crack tip enrichment value
            a3gp = sqrt(r)*sin(theta)*sin(theta/2);             % Node 3 crack tip enrichment value
            a4gp = sqrt(r)*sin(theta)*cos(theta/2);             % Node 4 crack tip enrichment value
            
            a1 = (a1gp - NODES(elementNodes(iN),5));                % Shifted alpha 1 enrichment value
            a2 = (a2gp - NODES(elementNodes(iN),7));                % Shifted alpha 2 enrichment value
            a3 = (a3gp - NODES(elementNodes(iN),9));                % Shifted alpha 3 enrichment value
            a4 = (a4gp - NODES(elementNodes(iN),11));               % Shifted alpha 4 enrichment value

            % Derivative of crack tip enrichment functions with respect to X
               Px = c*[-sin(theta/2)*ct             + cos(theta/2)*-st;...
                       cos(theta/2)*ct              + sin(theta/2)*-st;...
                      -sin(3*theta/2)*sin(theta)*ct + (sin(theta/2)+sin(3*theta/2)*cos(theta))*-st;...
                      -cos(3*theta/2)*sin(theta)*ct + (cos(theta/2)+cos(3*theta/2)*cos(theta))*-st];

               
               % Derivative of crack tip enrichment functions with respect to Y
               Py = c*[-sin(theta/2)*st             + cos(theta/2)*ct;...
                       cos(theta/2)*st              + sin(theta/2)*ct;...
                      -sin(3*theta/2)*sin(theta)*st + (sin(theta/2)+sin(3*theta/2)*cos(theta))*ct;...
                      -cos(3*theta/2)*sin(theta)*st + (cos(theta/2)+cos(3*theta/2)*cos(theta))*ct];

              B1 = [dNGlobal(1,iN)*a1+N(iN)*Px(1)          0;
                                       0            dNGlobal(2,iN)*a1+N(iN)*Py(1);
                              dNGlobal(2,iN)*a1+N(iN)*Py(1) dNGlobal(1,iN)*a1+N(iN)*Px(1)];

               B2 = [dNGlobal(1,iN)*a2+N(iN)*Px(2)          0;
                                       0            dNGlobal(2,iN)*a2+N(iN)*Py(2);
                              dNGlobal(2,iN)*a2+N(iN)*Py(2) dNGlobal(1,iN)*a2+N(iN)*Px(2)];

              B3 = [dNGlobal(1,iN)*a3+N(iN)*Px(3)          0;
                                       0            dNGlobal(2,iN)*a3+N(iN)*Py(3);
                              dNGlobal(2,iN)*a3+N(iN)*Py(3) dNGlobal(1,iN)*a3+N(iN)*Px(3)];

              B4 = [dNGlobal(1,iN)*a4+N(iN)*Px(4)          0;
                                       0            dNGlobal(2,iN)*a4+N(iN)*Py(4);
                              dNGlobal(2,iN)*a4+N(iN)*Py(4) dNGlobal(1,iN)*a4+N(iN)*Px(4)];

                        Bb = [B1 B2 B3 B4];
                        Benr(:,iB:(iB+7)) = Bb;
                        iB = iB+8;

                        if iGP == 1
                            index = 3;
                            for iAlpha = 1:4
                                Uenr(iLoc:(iLoc+1)) = [DISPLACEMENT(2*NODES(elementNodes(iN),iAlpha+index)-1,1) DISPLACEMENT(2*NODES(elementNodes(iN),iAlpha+index),1)];
                                index = index+1;
                                iLoc  = iLoc+2;
                            end
                        end
                   end
                 end

            stress = D*[Bu Benr]*[U;Uenr'];
            stressXX(iElem,iGP) = stress(1);
            stressYY(iElem,iGP) = stress(2);
            stressXY(iElem,iGP) = stress(3);
            stressVM(iElem,iGP) = sqrt(stress(1)^2+stress(2)^2-stress(1)*stress(2)+3*stress(3)^2);
        end
    end
end

% elemSxx = mean(stressXX, 2);
% elemSyy = mean(stressYY, 2);
% elemSxy = mean(stressXY, 2);
% elemSvm = mean(stressVM, 2);
% 
% % Option 1: Return element-based stresses (recommended for XFEM)
% % When plotting, use these element-average values.
% Sxx = elemSxx + Sxxr;
% Syy = elemSyy + Syyr;
% Sxy = elemSxy + Sxyr;
% Svm = elemSvm;

% Average nodal stress values
Sxx = zeros(nNodes,2); Syy = Sxx; Sxy = Sxx; Svm = Sxx;

    % Construct stress vectors
    for iE = 1:(nElements)
        for iN = 1:4
            nNodes = connectivity(iE,iN);
            Sxx(nNodes,:) = Sxx(nNodes,:) + [stressXX(iE,iN) 1];
            Sxy(nNodes,:) = Sxy(nNodes,:) + [stressXY(iE,iN) 1];
            Syy(nNodes,:) = Syy(nNodes,:) + [stressYY(iE,iN) 1];
            Svm(nNodes,:) = Svm(nNodes,:) + [stressVM(iE,iN) 1];
        end
    end

    % Average nodal stress values
    Sxx(:,1) = Sxx(:,1)./Sxx(:,2); Sxx(:,2) = [];
    Sxy(:,1) = Sxy(:,1)./Sxy(:,2); Sxy(:,2) = [];
    Syy(:,1) = Syy(:,1)./Syy(:,2); Syy(:,2) = [];
    Svm(:,1) = Svm(:,1)./Svm(:,2); Svm(:,2) = [];

    Sxx=Sxx+Sxxr;
    Syy=Syy+Syyr;
    Sxy=Sxy+Sxyr;

S1 = zeros(nNodes,1);
S2 = zeros(nNodes,1);
theta_p = zeros(nNodes,1);  % principal direction [rad]

for i = 1:nNodes
    sxx = Sxx(i);
    syy = Syy(i);
    sxy = Sxy(i);

    % Average and radius
    avg = 0.5 * (sxx + syy);
    radius = sqrt(((sxx - syy)/2)^2 + sxy^2 );

    % Principal stresses
    sigma1 = avg + radius;
    sigma2 = avg - radius;

    S1(i,1) = min(sigma1,sigma2);
    S2(i,1) = max(sigma1,sigma2);

    % Principal direction (radians)
    theta_p(i) = 0.5 * atan2(2*sxy, sxx - syy);  % angle from x-axis
end

% Optional: convert angle to degrees if needed
PP = theta_p * (180/pi);
                                                             
end