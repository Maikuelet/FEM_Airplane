function [Kel,KG,R] = StiffnessMatrix(El_L,T2,mat)

InputData;


%% Problem Dimensions

Ndim = size(x,1)*2;           % DoF for each node
Nnodes = size(x,2);           % Number of nodes
Nelements = size(Tnod,2);     % Number of elements
NnodesXelement = size(Tnod,1);% Number of nodes per element
Ndofs = Ndim*Nnodes;          % Total DoF of the system


%% Element Stiffness matrices

Kel = zeros(NnodesXelement*Ndim, NnodesXelement*Ndim, Nelements);
R = zeros(NnodesXelement*Ndim, NnodesXelement*Ndim, Nelements);
for e = 1:Nelements
    
    % ELEMENT ROTATION MATRIX
    
    alp=dat(2,e);   % Alpha 
    bet=dat(3,e);   % Beta
    gam=dat(4,e);   % Gamma
   
    R(:,:,e)=[cos(bet)*cos(gam),                              cos(bet)*sin(gam),                           sin(bet), 0,  0,  0,  0,  0,  0,  0,  0,  0;
       -(sin(alp)*sin(bet)*cos(gam))-(cos(alp)*sin(gam)),   -(sin(alp)*sin(bet)*sin(gam))+(cos(alp)*cos(gam)),  sin(alp)*cos(bet), 0,  0,  0,  0,  0,  0,  0,  0,  0;
       -(cos(alp)*sin(bet)*cos(gam))+(sin(alp)*sin(gam)),   -(cos(alp)*sin(bet)*sin(gam))-(sin(alp)*cos(gam)),  cos(alp)*cos(bet),0,  0,  0,  0,  0,  0,  0,  0,  0;
       0,  0,  0,cos(bet)*cos(gam),                              cos(bet)*sin(gam),                           sin(bet),  0,  0,  0,  0,  0,  0,;
       0,  0,  0,-(sin(alp)*sin(bet)*cos(gam))-(cos(alp)*sin(gam)),   -(sin(alp)*sin(bet)*sin(gam))+(cos(alp)*cos(gam)),  sin(alp)*cos(bet),  0,  0,  0,  0,  0,  0,;
       0,  0,  0,-(cos(alp)*sin(bet)*cos(gam))+(sin(alp)*sin(gam)),   -(cos(alp)*sin(bet)*sin(gam))-(sin(alp)*cos(gam)),  cos(alp)*cos(bet),  0,  0,  0,  0,  0,  0,;
       0,  0,  0,  0,  0,  0,   cos(bet)*cos(gam),                              cos(bet)*sin(gam),                           sin(bet),  0,  0,  0,;
       0,  0,  0,  0,  0,  0,   -(sin(alp)*sin(bet)*cos(gam))-(cos(alp)*sin(gam)),   -(sin(alp)*sin(bet)*sin(gam))+(cos(alp)*cos(gam)),  sin(alp)*cos(bet),  0,  0,  0,;
       0,  0,  0,  0,  0,  0,   -(cos(alp)*sin(bet)*cos(gam))+(sin(alp)*sin(gam)),   -(cos(alp)*sin(bet)*sin(gam))-(sin(alp)*cos(gam)),  cos(alp)*cos(bet),  0,  0,  0,;
       0,  0,  0,  0,  0,  0,  0,  0,  0,   cos(bet)*cos(gam),                              cos(bet)*sin(gam),                           sin(bet);
       0,  0,  0,  0,  0,  0,  0,  0,  0,   -(sin(alp)*sin(bet)*cos(gam))-(cos(alp)*sin(gam)),   -(sin(alp)*sin(bet)*sin(gam))+(cos(alp)*cos(gam)),  sin(alp)*cos(bet);
       0,  0,  0,  0,  0,  0,  0,  0,  0,   -(cos(alp)*sin(bet)*cos(gam))+(sin(alp)*sin(gam)),   -(cos(alp)*sin(bet)*sin(gam))-(sin(alp)*cos(gam)),  cos(alp)*cos(bet);
       ]; 
    
   % LOCAL ELEMENT SIFFNESS MATRIX 
   
   [Ae, Iye, Ize, Je] = BeamParameters(mat,dat,Tmat, e);
   
   le = El_L(e);
   Ee = mat(1,Tmat(e));
   Ge = mat(4,Tmat(e));
   
   Ke =[((Ee*Ae)/le),    0,      0,      0,      0,      0,      (-(Ee*Ae)/le),    0,      0,      0,      0,      0;
            0, ((12*Ee*Ize)/(le^3)),    0,      0,      0,  ((6*Ee*Ize)/(le^2)),    0, ((-12*Ee*Ize)/(le^3)),    0,      0,      0,  ((6*Ee*Ize)/(le^2));
            0, 0,   ((12*Ee*Iye)/(le^3)),    0,  ((-6*Ee*Iye)/(le^2)),    0, 0, 0, ((-12*Ee*Iye)/(le^3)),    0,  ((-6*Ee*Iye)/(le^2)),    0;
            0,  0,  0, ((Ge*Je)/le),   0,      0,      0,      0,      0,  ((-Ge*Je)/le),  0,      0;
            0, 0,   ((-6*Ee*Iye)/(le^2)),    0,  ((4*Ee*Iye)/(le)),    0, 0, 0, ((6*Ee*Iye)/(le^2)),    0,  ((2*Ee*Iye)/(le)),    0;
            0, ((6*Ee*Ize)/(le^2)),    0,      0,      0,  ((4*Ee*Ize)/(le)),    0, ((-6*Ee*Ize)/(le^2)),    0,      0,      0,  ((2*Ee*Ize)/(le));
            ((-Ee*Ae)/le),    0,      0,      0,      0,      0,      ((Ee*Ae)/le),    0,      0,      0,      0,      0;
            0, ((-12*Ee*Ize)/(le^3)),    0,      0,      0,  ((-6*Ee*Ize)/(le^2)),    0, ((12*Ee*Ize)/(le^3)),    0,      0,      0,  ((-6*Ee*Ize)/(le^2));
            0, 0,   ((-12*Ee*Iye)/(le^3)),    0,  ((6*Ee*Iye)/(le^2)),    0, 0, 0, ((12*Ee*Iye)/(le^3)),    0,  ((6*Ee*Iye)/(le^2)),    0;
            0,  0,  0, ((-Ge*Je)/le),   0,      0,      0,      0,      0,  ((Ge*Je)/le),  0,      0;
            0, 0,   ((-6*Ee*Iye)/(le^2)),    0,  ((2*Ee*Iye)/(le)),    0, 0, 0, ((6*Ee*Iye)/(le^2)),    0,  ((4*Ee*Iye)/(le)),    0;
            0, ((6*Ee*Ize)/(le^2)),    0,      0,      0,  ((2*Ee*Ize)/(le)),    0, ((-6*Ee*Ize)/(le^2)),    0,      0,      0,  ((4*Ee*Ize)/(le));
            ];
    
   % GLOBAL ELEMENT STIFFNESS MATRIX
   Ke = R(:,:,e).' * Ke * R(:,:,e);
   
   % ELEMENT MATRIX ASSEMBLY
   for r=1:(NnodesXelement*Ndim)
        for s=1:(NnodesXelement*Ndim)
            Kel(r,s,e)=Ke(r,s);
        end
   end
    
end

%% GLOBAL STIFFNESS MATRIX
KG=zeros(Ndofs,Ndofs);
for e=1:Nelements
    for i=1:(NnodesXelement*Ndim)
        I = T2(i,e);
        for j=1:(NnodesXelement*Ndim)
            J = T2(j,e);
            KG(I,J)=KG(I,J)+ Kel(i,j,e);
        end
    end
end


end