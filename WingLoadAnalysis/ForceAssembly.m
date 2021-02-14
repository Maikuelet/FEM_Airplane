function [Fext,Drag_T,Lift_T,El_W] = ForceAssembly(El_L,rho_eff,Lift,Drag,R,mat,T2)
    
InputData;


%% Problem Dimensions

Ndim = size(x,1)*2;           % DoF for each node
Nnodes = size(x,2);           % Number of nodes
Nelements = size(Tnod,2);     % Number of elements
NnodesXelement = size(Tnod,1);% Number of nodes per element
Ndofs = Ndim*Nnodes;          % Total DoF of the system



%% FORCE VECTOR AND MATRIX ASSEMBLY
El_W = zeros(Nelements,1);  % Weight Vector
El_T = zeros(Nelements,1);  % Tension Vector
fe_local = zeros(Ndim*NnodesXelement,1);    % External force local
fe = zeros(Ndim*NnodesXelement,Nelements);  % External force global element
Fext=zeros(Ndofs,1);                        % System External force vector

Lift_T = 0;
Drag_T = 0;


for e = 1:Nelements
    Drag_T = Drag_T + Drag(e)*El_L(e);
    Lift_T = Lift_T + Lift(e)*El_L(e);
end

for e = 1:Nelements
 
    le = El_L(e);
    [Ae] = BeamParameters(mat,dat,Tmat, e);
    
    % Element distributed Weight
    El_W(e) = rho_eff(Tmat(e))*Ae*g;
    
    if e==7 || e==8
      El_W(e) = El_W(e) + ((M_eng*g)/(2*le));
      El_T(e) = Drag_T / (2*le);
    end
    
    qx = (Drag(e)- El_T(e));
    qz = (Lift(e)-El_W(e));
    qe = [qx; 0; qz; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
    
    % COMPONENTS IN LOCAL COORDINATE SYSTEM
    qe_local = R(:,:,e)*qe;

%     EQUIVALENT ELEMENT FORCE VECTOR IN LOCAL COORDINATES
    for d = 1:Ndim*NnodesXelement
        % Axial distributed load x-direction
        if d == 1 || d==7
            ind = 1;
            fe_local(d) = (qe_local(ind)*le/2);
        end
        % Axial distributed load y-direction
        if d == 2 || d==6 || d==8 || d==12
            ind = 2;
            if d == 2 || d==8
                fe_local(d) = (qe_local(ind)*le/2)*1;
            end
            if d == 6
                fe_local(d) = (qe_local(ind)*le/2)* (le/6);
            end
            if d == 12
                fe_local(d) = (qe_local(ind)*le/2)* (-le/6);
            end
        end
        % Axial ditributed load z-direction
        if d == 3 || d==5 || d==9 || d==11
            ind = 3;
            if d == 3 || d==9
                fe_local(d) = (qe_local(ind)*le/2)*1;
            end
            if d == 5
                fe_local(d) = (qe_local(ind)*le/2)*(-le/6);
            end
            if d == 11
                fe_local(d) = (qe_local(ind)*le/2)*(le/6);
            end
        end
    end
    
    % EQUIVALENT ELEMENT FORCE VECTOR IN GLOBAL COORDINATES
    fe(:,e) = R(:,:,e).' * fe_local;
    
    % FORCE ASSEMBLY
    for i = 1:NnodesXelement*Ndim
        I = T2(i,e);
        Fext(I)=Fext(I)+fe(i,e);
    end   
end


end