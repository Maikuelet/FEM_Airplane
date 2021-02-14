function [El_L,El_M,M_spar,M_rib, V_spar, V_rib,rho_eff,M_struct,V_tot] = StructureMass(mat)

InputData;


%% Problem dimensions

Ndim = size(x,1)*2;           % DoF for each node
Nnodes = size(x,2);           % Number of nodes
Nelements = size(Tnod,2);     % Number of elements
NnodesXelement = size(Tnod,1);% Number of nodes per element
Ndofs = Ndim*Nnodes;          % Total DoF of the system

%% Structure Computations

El_L = zeros (Nelements,1);     %Save Element Length
El_M = zeros (Nelements,1);     %Save Element Mass


M_spar = 0;
M_rib  = 0;
V_rib  = 0;
V_spar = 0;
V_tot  = 0;


for e = 1:Nelements
    x1=x(1,Tnod(1,e));
    y1=x(2,Tnod(1,e));
    z1=x(3,Tnod(1,e));
    x2=x(1,Tnod(2,e));
    y2=x(2,Tnod(2,e));
    z2=x(3,Tnod(2,e));
    
    El_L(e) = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2);
    
    [Ae] = BeamParameters(mat,dat,Tmat, e);
    
    El_V  = Ae*El_L(e);
    rho_e = mat(3,Tmat(e));
    El_M(e) =El_V*rho_e;
    
    if Tmat(e)==1
       M_spar = El_M(e) + M_spar;
       V_spar = El_V + V_spar;
    end
    if Tmat(e)== 2 % Rib
       M_rib = El_M(e) + M_rib;
       V_rib = El_V + V_rib;
    end  
    
    V_tot = El_V + V_tot;

end

rho_hat = (M_w - M_spar-M_rib)/ (V_rib+V_spar);
rho_eff(1)  =  mat(3,1) + rho_hat; % Spar
rho_eff(2)  =  mat(3,2) + rho_hat; % Rib

M_struct = sum(El_M);
% Mtot  = V_rib*rho_eff_rib + V_spar*rho_eff_spar + M_w;
% Mtot2 = M_spar + M_rib + M_w;
%disp('Total Mass 1'); disp(Mtot);
%disp('Total Mass 2'); disp(Mtot2);



end