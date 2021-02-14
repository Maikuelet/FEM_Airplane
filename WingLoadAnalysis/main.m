%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%    FINITE ELEMENT METHOD     %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     STRUCTURAL ANALYSIS      %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  Miquel Altadill Llasat  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code developed for the structural analysis of a commercial 
% aircraft wing.  

clc
clear all
close all

%% =========================== INPUT DATA ============================== %%
%InputData
% - Inside each function (Cleaner Workspace)

%% ==========================  PREPROCESS =============================  %%
[T2,Fnod, mat,Lift,Drag] = Preprocess();
[El_L,El_M,M_spar,M_rib, V_spar, V_rib,rho_eff,M_struct,V_tot] = StructureMass(mat);

%% ============================= SOLVER ===============================  %%
[Kel,KG,R] = StiffnessMatrix(El_L,T2,mat);
[Fext,Drag_T,Lift_T,El_W] = ForceAssembly(El_L,rho_eff,Lift,Drag,R,mat,T2);
[U,F_int,u_int,N_x,Q_y,Q_z,T_x,M_y,M_z,rf,Rr] = solver(Fnod,KG,Kel,R,T2,Fext);


%% ========================== POSTPROCESS =============================  %%
PostProcess(M_struct,M_rib,M_spar,V_rib,V_spar,Lift_T,Drag_T,V_tot,U,u_int,...
    F_int,El_L,N_x,Q_y,Q_z,T_x,M_y,M_z,R,Rr,El_W);
