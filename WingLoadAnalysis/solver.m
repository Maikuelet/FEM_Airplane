function [U,F_int,u_int,N_x,Q_y,Q_z,T_x,M_y,M_z,rf,Rr] = solver(Fnod,KG,Kel,R,T2,Fext)

InputData;

%% Problem Dimensions

Ndim = size(x,1)*2;           % DoF for each node
Nnodes = size(x,2);           % Number of nodes
Nelements = size(Tnod,2);     % Number of elements
NnodesXelement = size(Tnod,1);% Number of nodes per element
Ndofs = Ndim*Nnodes;          % Total DoF of the system



%% SYSTEM DEGREES OF FREEDOM

% RESTRICTED DoF's:
    
NFnod = size (Fnod,2);
vR = zeros(1,NFnod);    % Restricted DoF global
uR = zeros(NFnod,1);    % Restricted displacements

for i = 1:NFnod
   uR(i,1) = Fnod(3,i); 
end
for i = 1:NFnod
   node = Fnod(1,i); 
   dof  = Fnod(2,i);
   vR(1,i) = Ndim*node - Ndim + dof;
end


% FREE DoF's
vL = setdiff(1:Ndofs,vR);



%% REDUCED SYSTEM SOLUTION

KLL = KG(vL,vL);
KLR = KG(vL,vR);
KRL = KG(vR,vL);
KRR = KG(vR,vR);
fextL = Fext(vL,1);
fextR = Fext(vR,1);

% LINEAR SYSTEM SOLUTION
ul = KLL\(fextL - KLR*uR);
Rr = KRR*uR + KRL*ul - fextR;

% GLOBAL DISPLACEMENTS ASSEMBLY
U = zeros(Ndofs,1);
U(vL) = ul;
U(vR) = uR;

% GLOBAL REACTIONS ASSEMBLY
rf = zeros(Ndofs,1);
rf(vL) = 0;
rf(vR) = Rr;


%% INTERNAL FORCES AND MOMENTS

u_int = zeros(NnodesXelement*Ndim, Nelements);  % Internal displacements and rot
F_int = zeros(NnodesXelement*Ndim, Nelements);  % Internal forces in local
N_x = zeros(1,Nelements);
Q_y = zeros(1,Nelements);
Q_z = zeros(1,Nelements);
T_x = zeros(1,Nelements);
M_y = zeros(2,Nelements);
M_z = zeros(2,Nelements);

for  e = 1:Nelements
    for i = 1:NnodesXelement*Ndim
        I=T2(i,e);
        % INTERNAL GLOBAL DISPLACEMENTS AND ROTATIONS 
        u_e(i,1) = U(I,1);
    end
    % INTERNAL LOCAL DISPLACEMENTS AND ROTATIONS 
    u_int(:,e) = R(:,:,e) * u_e;

    % INTERNAL FORCE VECTOR IN GLOBAL
    F_e = Kel(:,:,e) * u_e;
    % INTERNAL FORCE VECTOR IN LOCAL
    F_int(:,e) = R(:,:,e) * F_e;
    % AXIAL FORCE
    N_x(e) = -F_int(1,e);  % X direction
    % SHEAR FORCE 
    Q_y(e) = - F_int(2,e);     % Y direction
    Q_z(e) = - F_int(3,e);     % Z direction
    %TORSION MOMENT
    T_x(e) = - F_int(4,e);      % Arround X-axis
    % BENDING MOMENT (Linear)
    M_y(1,e)  = -F_int(5,e); 
    M_y(2,e)  = F_int(11,e);     % Arround Y-axis
    M_z(1,e)  = -F_int(6,e);
    M_z(2,e)  =  F_int(12,e);     % Arround Z-axis
    
end


end