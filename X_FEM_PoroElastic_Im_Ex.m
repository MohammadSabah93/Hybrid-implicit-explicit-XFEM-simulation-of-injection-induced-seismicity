clc;
clear;
close all;


%% 2D Mesh Generation 

% Generate 2D mesh
domain_length = 100;
domain_width = 100;
nex = 100;
ney = 100;

crackStart = [20 ,20.5];
crackEnd = [80, 80.5];

CRACK=[crackStart;crackEnd];

[connectivity, node_coordinates, nNodes, nElements] = MeshGeneration_2D(domain_length, domain_width, nex, ney, crackStart, crackEnd);

[omega, PSI, PHI, NODES, interface_elements, interface_nodes]=levelSet(connectivity, node_coordinates, nNodes, domain_length, nex, ney, CRACK);

nodesPerElement = size(connectivity, 2);
if nodesPerElement == 4
   elementType = 'QUAD4';
else
   elementType = 'TRI3';
end

%% Calculate unit normal and tangential vectors on fault

  % Extract coordinates of the crack
    x1 = crackStart(1);
    y1 = crackStart(2);
    x2 = crackEnd(1);
    y2 = crackEnd(2);

    % Compute the direction vector and length of the crack
    crackVector = [x2 - x1; y2 - y1];
    crackLengthSquared = dot(crackVector, crackVector); % Avoid recomputing square

    % Normalize the normal vector to the crack
    n = [crackVector(2); -crackVector(1)] / sqrt(crackLengthSquared);

    % Normalize the tangential vector to the crack
    t = crackVector / sqrt(crackLengthSquared);

%% Model Parameters

params = defineModelParameters();

%% Boundary Conditions

boundaryConditions = defineBoundaryConditions(node_coordinates,connectivity);

%% Pre-allocate solution vectors
[DOF_u,DOF_p,DISP] = calcDOF(NODES);
nInterfaceNodes = max(interface_elements(:));
nInterfaceElem = size(interface_elements,1);
a = zeros(2*nnz(NODES(:,2)),1);
Ux = zeros(nNodes, 1); 
Uy = zeros(nNodes, 1); 
U_ddotx = zeros(nNodes, 1); 
U_ddoty = zeros(nNodes, 1);
U_ddot_total = zeros(nNodes, 1);
U_dotx = zeros(nNodes, 1); 
U_doty = zeros(nNodes, 1);
U_dot_total = zeros(nNodes, 1);
P = zeros(DOF_p, 1);
Pp = zeros(nNodes, 1);
Dn = zeros(nInterfaceElem,1); 
Ds = zeros(nInterfaceElem,1); 
Sn = zeros(nInterfaceElem,1); 
Ss = zeros(nInterfaceElem,1);
Pn = zeros(nInterfaceElem,1); 
W = zeros(nInterfaceElem,1);
LM = zeros(nInterfaceNodes,1);
slip = zeros(nInterfaceElem,1);
uf = zeros(nInterfaceElem,1);
theta = zeros(nInterfaceElem,1);
Vs = zeros(nInterfaceElem,1);
stickflag = zeros(nInterfaceElem,1);
duf = zeros(nInterfaceElem,1);

%% Initial conditions

Nt = 10000;
dt_c=3600;
dt = dt_c;
W (:,1) = params.W0;
uf(:,1) = params.uf0;
theta(:,1) = params.theta0;
P(:,1)=params.P0;
Time(1,1)=0;
dampingFactor = 1;
Ss(:,1)=-9.5e+6;
Sn(:,1)=-3e+7;
Ds(:,1)=-0.0009;
landa=0.6;
beta=0.5;

%% Assembling Process for PoroElasticity
globalM = MassMatrix(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, elementType, omega, domain_length, nex);
globalC = InertialMatrix(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, elementType, omega, domain_length, nex);

globalK = StiffnessMatrix(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, omega, domain_length, nex, elementType);
globalH = ConductanceMatrix(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, elementType);
globalS = StorageMatrix(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, elementType);
globalQ = CouplingMatrix(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, elementType,omega, domain_length, nex);
globalDa = dampingBoundary(boundaryConditions, node_coordinates, connectivity, NODES, elementType, params);


globalQ_inter = CouplingInterface(connectivity, node_coordinates, PSI, NODES, CRACK, n, elementType, PHI, omega, domain_length, nex);
globalFu_ext1 = TractionForces(boundaryConditions, node_coordinates, connectivity, NODES, elementType);
globalFu_ext2 = BodyForceSolid(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, omega, domain_length, nex, elementType);
globalFu_ext3 = InsituStress(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, omega, domain_length, nex, elementType);
globalFp_ext1 = BoundaryFluidSource(boundaryConditions, node_coordinates, connectivity, NODES, elementType);
globalFp_ext2 = BodyForceFluid(params, node_coordinates, connectivity, PSI, PHI, NODES, CRACK, elementType);
globalFp_ext3 = InjectionSource(boundaryConditions, NODES);
globalG_inter = Lagrange(connectivity, node_coordinates, NODES, CRACK, n, interface_elements, elementType, PHI, omega, domain_length, nex);
globalG_inter2 = StabLagrange(params,interface_elements,interface_nodes);
globalQ = globalQ + globalQ_inter;
globalFu = globalFu_ext1 + globalFu_ext2 - globalFu_ext3;

omega1=pi*(2*1-1)*params.cs/(2*domain_width);
omega3=pi*(2*3-1)*params.cs/(2*domain_width);
a_ray=2*omega1*omega3*0.04/(omega1+omega3);
b_ray=2*0.04/(omega1+omega3);
globalD = a_ray*globalM+b_ray*globalK + globalDa;

U = zeros(DOF_u, 1); 
U_ddot = zeros(DOF_u, 1); 
U_dot = zeros(DOF_u, 1); 
ResidualHistory = cell(Nt+1,1);

%% Apply Boundary conditions

fixedDOF_disp = boundaryConditions.fixedDOF_disp;
dispDOFs = fixedDOF_disp(:,1);
fixedPressureDOF = boundaryConditions.fixedPressureDOF;

if ~isempty(fixedPressureDOF)
    pressureDOFs = fixedPressureDOF(:,1) + DOF_u;
    fixedPressures = fixedPressureDOF(:,2);
else
    pressureDOFs = [];
    fixedPressures = [];
    diagPressIdx = [];
end

%% Newton-Raphson iteration

numElems = size(interface_elements, 1);
element_length = zeros(numElems, 1);  % [m]

    % Compute element lengths
    for e = 1:numElems
        n1 = interface_elements(e,1);
        n2 = interface_elements(e,2);
        coord1 = interface_nodes(n1, :);
        coord2 = interface_nodes(n2, :);
        element_length(e) = norm(coord2 - coord1);
    end

tolerance        = 1e-6;
maxIter          = 15;
useImplicit      = true;      % start in implicit mode
Vth              = 1e-3;      % <<< slip-rate threshold for switching
wasImplicit = useImplicit;

for nt = 1:20000

    % ==================================================================
    %   A) IMPLICIT MODE (current Newmark + Newton–Raphson)
    % ==================================================================
if useImplicit

    % Initial guess for next time step
    U(:,nt+1)  = U(:,nt);
    P(:,nt+1)  = P(:,nt);
    W(:,nt+1)  = W(:,nt);
    Sn(:,nt+1) = Sn(:,nt);
    Ss(:,nt+1) = Ss(:,nt);
    uf(:,nt+1) = uf(:,nt);
    stickflag(:,nt+1) = stickflag(:,nt);
    Vs(:,nt+1) = Vs(:,nt);

        if all(Vs(:,nt) == 1e-9)
            dt_ev = inf;
        else
            dt_ev = min(0.1*params.Dc ./ Vs(:,nt));
        end
        dt_candidate = min([max(1e-3, dt_ev), dt_c]);
        dt = min(dt_candidate,1.2*dt);

%% Apply Boundary conditions

fixedDOF_disp = boundaryConditions.fixedDOF_disp;
dispDOFs = fixedDOF_disp(:,1);
fixedPressureDOF = boundaryConditions.fixedPressureDOF;

if ~isempty(fixedPressureDOF)
    pressureDOFs = fixedPressureDOF(:,1) + DOF_u;
    fixedPressures = fixedPressureDOF(:,2);
else
    pressureDOFs = [];
    fixedPressures = [];
    diagPressIdx = [];
end
    iter = 0;
    Norm = Inf;

while Norm > tolerance && iter < maxIter

        iter = iter + 1;

        dK_dU = StiffnessInterface_Penalty(params, connectivity, node_coordinates, NODES, CRACK, n, t, Sn(:,nt+1), duf, dt, stickflag(:,nt+1), elementType, PHI, omega, domain_length, nex);
        globalF_inter = FaultForce(connectivity, node_coordinates, NODES, CRACK, Sn(:,nt+1), Ss(:,nt+1), t, n, elementType, PHI, omega, domain_length, nex);
        [H_inter, S_inter, Q_inter, F_inter] = interface_flow(params, connectivity, node_coordinates, PSI, NODES, CRACK, t, W(:,nt+1), elementType);
        H = globalH + H_inter;
        S = globalS + S_inter;

        % Form residuals 
        Gu = globalFu + globalM * (1/(beta*dt^2)*U(:,nt) + 1/(beta*dt)*U_dot(:,nt) + (1/(2*beta)-1)*U_ddot(:,nt)) + ...
                        globalD * (landa/(beta*dt) * U(:,nt) + (landa/beta-1) * U_dot(:,nt) + dt*(landa/(2*beta)-1)*U_ddot(:,nt));
        Gp = (globalFp_ext1+globalFp_ext2+globalFp_ext3) + F_inter + globalC * (1/(beta*dt^2)*U(:,nt) + 1/(beta*dt)*U_dot(:,nt) + (1/(2*beta)-1)*U_ddot(:,nt)) + ...
                        (globalQ'+Q_inter)* (landa/(beta*dt)*U(:,nt) + (landa/beta-1)*U_dot(:,nt) + dt * (landa/(2*beta)-1)*U_ddot(:,nt)) + S/dt * P(:,nt);
                      

        Residual  = [-landa/(beta*dt) * ((1/(beta*dt^2) * globalM + landa/(beta*dt) * globalD + globalK) * U(:,nt+1) - globalQ * P(:,nt+1) + globalF_inter - Gu);...
                    (1/(beta*dt^2) * globalC + landa/(beta*dt) * (globalQ'+Q_inter)) * U(:,nt+1) + (S/dt + H) * P(:,nt+1) - Gp];


        % Form jacobian
        K = globalK;
        K(1:DOF_u, 1:DOF_u) = K(1:DOF_u, 1:DOF_u) + dK_dU;
        Jacobian = [-landa/(beta*dt) * ((1/(beta*dt^2) * globalM + landa/(beta*dt) * globalD + K)), landa/(beta*dt) * globalQ; ...
                    (1/(beta*dt^2) * globalC + landa/(beta*dt) * (globalQ'+Q_inter)), (1/dt*S + H)];


        % Apply essential boundary conditions (displacement)
        Jacobian(dispDOFs, :) = 0;
        Jacobian(:, dispDOFs) = 0;
        diagDispIdx = sub2ind(size(Jacobian), dispDOFs, dispDOFs);
        Jacobian(diagDispIdx) = 1;
        Residual(dispDOFs) = 0;

        % Apply essential boundary conditions (pressure)
        if ~isempty(pressureDOFs)
        Jacobian(pressureDOFs, :) = 0;
        Jacobian(:, pressureDOFs) = 0;
        diagPressIdx = sub2ind(size(Jacobian), pressureDOFs, pressureDOFs);
        Jacobian(diagPressIdx) = 1;
        Residual(pressureDOFs) = P(fixedPressureDOF(:,1), nt+1) - fixedPressures;
        end

% Apply scaling to J and r
Jt = sparse(Jacobian);  
rt = -Residual;
deltaX = Jt \ rt;

% Normalize Displacement and Pressure Corrections**
        ResidualU = norm(Residual(1:(DOF_u),1)) / norm(globalFu);
        ResidualP = norm(Residual(DOF_u+1:end,1)) / (norm(globalFp_ext3)+1e-5);
% ResidualU = norm(Residual(1:(DOF_u),1));
%         ResidualP = norm(Residual(DOF_u+1:end,1));
        % Compute L2 Norm with Normalized Values**
        Norm = sqrt(norm(ResidualU, 2)^2 + norm(ResidualP, 2)^2);
        % ResidualHistory{nt+1,1}=[ResidualHistory{nt+1}; Norm];
        fprintf('Time step %d, Iteration %d: Norm = %e\n', nt, iter, Norm);

        % Update displacement and pressure for time step nt+1
        U(:, nt+1) = U(:, nt+1) + dampingFactor*deltaX(1:DOF_u);  % update displacement
        P(:, nt+1) = P(:, nt+1) + dampingFactor*deltaX(DOF_u+1:end);  % Update pressure

        Ux(:,nt+1) = U(1:2:2*nNodes,nt+1);
        Uy(:,nt+1) = U(2:2:2*nNodes,nt+1);
        a(:,nt+1) = U(2*nNodes+1:2*nNodes+2*nnz(NODES(:,2)),nt+1);
        % LM(:,nt+1) = U(DOF_u+1:end,nt+1);
        Pp(:,nt+1) = P(1:nNodes,nt+1);

        [pn,pt,gn,gt,~] = Interface(params, connectivity, node_coordinates, NODES, CRACK, U(:,nt+1), n, t, interface_elements, LM(:,1), elementType, omega, domain_length, nex, PHI);
        Dn(:,nt+1) = gn;
        Ds(:,nt+1) = gt;
        Sn(:,nt+1) = pn;
        Ss(:,nt+1) = pt;
        W(:,nt+1) = params.W0 + Dn(:,nt+1);

        [Ss, slip, Vs, uf, duf, theta, stickflag] = updateInterface(params, Ds, Ss, Sn, Vs, theta, uf, slip, stickflag, nt, dt, iter);

        
        
        if Norm < tolerance
            fprintf('Convergence achieved at time step %d, iteration %d\n', nt, iter);
            break;
        end
end
U_ddot(1:DOF_u,nt+1) = 1/(beta*dt^2) * (U(1:DOF_u, nt+1)-U(1:DOF_u, nt)) - 1/(beta*dt) * U_dot(1:DOF_u,nt) - (1/(2*beta)-1) * U_ddot(1:DOF_u,nt);
U_ddotx(:,nt+1) = U_ddot(1:2:2*nNodes,nt+1);
U_ddoty(:,nt+1) = U_ddot(2:2:2*nNodes,nt+1);
U_ddot_total(:,nt+1) = sqrt(U_ddotx(:,nt+1).^2 + U_ddoty(:,nt+1).^2);
U_dot(1:DOF_u,nt+1) = landa/(beta*dt)*(U(1:DOF_u,nt+1)-U(1:DOF_u,nt)) - (landa/beta-1)*U_dot(1:DOF_u,nt) - dt * (landa/(2*beta)-1)*U_ddot(1:DOF_u,nt);
U_dotx(:,nt+1) = U_dot(1:2:2*nNodes,nt+1);
U_doty(:,nt+1) = U_dot(2:2:2*nNodes,nt+1);
U_dot_total(:,nt+1)  = sqrt(U_dotx(:,nt+1).^2  + U_doty(:,nt+1).^2);
% ==================================================================
    %   B) EXPLICIT MODE (after Vdyn exceeds threshold)
    % ==================================================================
else

    
    globalF_inter = FaultForce(connectivity, node_coordinates, NODES, CRACK, Sn(:,nt), Ss(:,nt), t, n, elementType, PHI, omega, domain_length, nex);

    % ---- Known half-step velocity u_dot_{i-1/2} from implicit result ----
    % % u_dot_{i-1/2} = u_dot_i - (dt/2)*u_ddot_i          (Eq. 17.3.8, backward)
    Udot_half_prev = U_dot(1:DOF_u, nt) ...
                      - 0.5 * dt * U_ddot(1:DOF_u, nt);

    % ---- Effective mechanical load f_i at time t_i on displacement DOFs ----
    F_u = globalFu(1:DOF_u) ...
          - globalF_inter ...
          + globalQ(1:DOF_u,:) * P(:,nt);

    % ---- Eq. of motion at time t_i (17.3.7) with substitution (17.3.8):
    % (M + dt/2 * C) * u_ddot_i = F_i - C*u_dot_{i-1/2} - K*u_i   (17.3.9)
    A   = globalM(1:DOF_u,1:DOF_u) ...
        + 0.5 * dt * globalD(1:DOF_u,1:DOF_u);

    rhs = F_u ...
          - globalD(1:DOF_u,1:DOF_u) * Udot_half_prev...
          - globalK(1:DOF_u,1:DOF_u) * U(1:DOF_u, nt);

    % ---- Apply displacement Dirichlet BCs on (M + dt/2 C) ----
    A(dispDOFs, :) = 0;
    A(:, dispDOFs) = 0;
    diagDispIdx    = sub2ind(size(A), dispDOFs, dispDOFs);
    A(diagDispIdx) = 1;
    rhs(dispDOFs)  = 0;

    % ---- Solve for acceleration u_ddot_i at time t_i ----
    U_ddot(1:DOF_u, nt) = A \ rhs;
    U_ddot(dispDOFs, nt) = 0;

    % ---- Update half-step velocity: u_dot_{i+1/2} = u_dot_{i-1/2} + dt*u_ddot_i  (17.3.6b)
    Udot_half_next = Udot_half_prev + dt * U_ddot(1:DOF_u, nt);

    % ---- Update displacement to next step:
    % u_{i+1} = u_i + dt * u_dot_{i+1/2}                                    (17.3.6a)
    U(1:DOF_u, nt+1) = U(1:DOF_u, nt) + dt * Udot_half_next;
    U(dispDOFs, nt+1) = 0;

    U_dot(1:DOF_u, nt+1)  = Udot_half_next;
    U_dot(dispDOFs, nt+1) = 0;

    U_ddot(1:DOF_u, nt+1) = U_ddot(1:DOF_u, nt);
    U_ddot(dispDOFs, nt+1) = 0;
    Udot_half_prev=Udot_half_next;
    
    [pn,pt,gn,gt,~] = Interface(params, connectivity, node_coordinates, NODES, CRACK, U(:,nt+1), n, t, interface_elements, LM(:,1), elementType, omega, domain_length, nex, PHI);
        Dn(:,nt+1) = gn;
        Ds(:,nt+1) = gt;
        Sn(:,nt+1) = pn;
        Ss(:,nt+1) = pt;
        W(:,nt+1) = params.W0 + Dn(:,nt+1);

    [Ss, slip, Vs, uf, duf, theta, stickflag] = updateInterface_expilict(params, Ds, Ss, Sn, Vs, theta, uf, slip, stickflag, nt, dt);


U_ddotx(:,nt+1) = U_ddot(1:2:2*nNodes,nt+1);
U_ddoty(:,nt+1) = U_ddot(2:2:2*nNodes,nt+1);
U_ddot_total(:,nt+1) = sqrt(U_ddotx(:,nt+1).^2 + U_ddoty(:,nt+1).^2);
U_dotx(:,nt+1) = U_dot(1:2:2*nNodes,nt+1);
U_doty(:,nt+1) = U_dot(2:2:2*nNodes,nt+1);
U_dot_total(:,nt+1)  = sqrt(U_dotx(:,nt+1).^2  + U_doty(:,nt+1).^2);

 % 2) FLUID STEP (Backward Euler in time using u_dot^{n+1}, u_ddot^{n+1})
 [H_inter, S_inter, Q_inter, F_inter] = interface_flow(params, connectivity, node_coordinates, PSI, NODES, CRACK, t, W(:,nt+1), elementType);
    H = globalH + H_inter;
    S = globalS + S_inter;
    %------------------------------------------------------------------
    Fp = (globalFp_ext1 + globalFp_ext2 + globalFp_ext3) ...
         - (globalQ' + Q_inter) * U_dot(:, nt+1) ...
         - globalC * U_ddot(1:DOF_u, nt+1);
    
% Backward Euler: (S + dt*H) P^{n+1} = S*P^{n} + dt*Fp^{n+1}
A_p = S + dt * H;
B_p = S * P(:,nt) + dt * Fp;

% --- Apply pressure Dirichlet BCs correctly ---
if ~isempty(pressureDOFs)
    A_p(pressureDOFs, :) = 0;
    A_p(:, pressureDOFs) = 0;
    diagPressIdx         = sub2ind(size(A_p), pressureDOFs, pressureDOFs);
    A_p(diagPressIdx)    = 1;

    % Enforce prescribed pressure directly
    B_p(pressureDOFs) = fixedPressures;
end

% --- Solve for P^{n+1} ---
P(:, nt+1) = A_p \ B_p;

        % ---- Post-process interface fields, Vs, etc. as usual ----
        Ux(:,nt+1) = U(1:2:2*nNodes,nt+1);
        Uy(:,nt+1) = U(2:2:2*nNodes,nt+1);
        a(:,nt+1)  = U(2*nNodes+1:2*nNodes+2*nnz(NODES(:,2)),nt+1);
        Pp(:,nt+1) = P(1:nNodes,nt+1);  
        fprintf('Explicit time step: %d\n', nt);

end

 % ==============================================================
    %   C) MODE SWITCHING LOGIC (IMPLICIT <-> EXPLICIT)
    % ==============================================================

    % Peak slip velocity at *current* step
    Vs_step = Vs(:, nt+1);
    Vdyn    = max(Vs_step);

    % Decide which mode to use for the *next* time step
    if Vdyn > Vth
        useImplicit_next = false;   % some element is dynamic → explicit
         
    else
        useImplicit_next = true;    % all quiet → implicit
    end

    % Detect transitions and adjust dt / flags
    if wasImplicit && ~useImplicit_next
        % IMPLICIT -> EXPLICIT
        fprintf('*** Switching to EXPLICIT at time step %d, Vdyn = %e m/s ***\n', nt+1, Vdyn);
        dt=1e-4;
        pressureDOFs = fixedPressureDOF(:,1);
        wasImplicit = false;
    elseif ~wasImplicit && useImplicit_next
        % EXPLICIT -> IMPLICIT
        fprintf('*** Switching back to IMPLICIT at time step %d, Vdyn = %e m/s ***\n', nt+1, Vdyn);
         U_dot(1:DOF_u, nt+1) = (U(1:DOF_u, nt+1) - U(1:DOF_u, nt)) / dt;
         globalF_inter = FaultForce(connectivity, node_coordinates, NODES, CRACK, Sn(:,nt+1), Ss(:,nt+1), t, n, elementType, PHI, omega, domain_length, nex);
         A   = globalM(1:DOF_u,1:DOF_u);
         F_u = globalFu(1:DOF_u) - globalF_inter + globalQ(1:DOF_u,:) * P(:,nt+1);
         rhs = F_u - globalD(1:DOF_u,1:DOF_u) * U_dot(1:DOF_u, nt+1) - globalK(1:DOF_u,1:DOF_u) * U(1:DOF_u, nt+1);
         A(dispDOFs, :) = 0;
         A(:, dispDOFs) = 0;
         diagDispIdx    = sub2ind(size(A), dispDOFs, dispDOFs);
         A(diagDispIdx) = 1;
         rhs(dispDOFs)  = 0;
         U_ddot(1:DOF_u, nt+1) = A \ rhs;
         U_ddot(dispDOFs, nt+1) = 0;

fixedDOF_disp = boundaryConditions.fixedDOF_disp;
dispDOFs = fixedDOF_disp(:,1);
fixedPressureDOF = boundaryConditions.fixedPressureDOF;

if ~isempty(fixedPressureDOF)
    pressureDOFs = fixedPressureDOF(:,1) + DOF_u;
    fixedPressures = fixedPressureDOF(:,2);
else
    pressureDOFs = [];
    fixedPressures = [];
    diagPressIdx = [];
end
    wasImplicit = true;

    end

    useImplicit = useImplicit_next; % update mode for next iteration


Time(nt+1,1) = dt + Time(nt,1);    
current_time_hours = Time(end) / 3600;
if current_time_hours >= 1000
   fprintf('Reached total simulation time of 1000 hours. Quitting simulation.\n');
   break;
end
end

%% Calculate seismicity parameters

[AS, M0, Mw, CFS, SSD, ST, seisEff] = SeismicityParameters(interface_elements, interface_nodes, Vs, slip, Ss, Sn, uf, params);

%% Calculate stresses 

[Sxx,Sxy,Syy,Svm,S1,S2,theta_p] = elemStress2(params, U, connectivity, node_coordinates, PSI, NODES, CRACK, omega, domain_length, nex, elementType, Nt);

%% Post-Processing and Visualization

postProcess(node_coordinates, U(1:2*nNodes,end), Pp(:,end), connectivity, Sxx(:,end), Syy(:,end), Sxy(:,end), elementType);