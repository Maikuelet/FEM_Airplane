function [T2,Fnod, mat,Lift,Drag] = Preprocess()

%% Input Data

InputData;

% Material properties storage
mat = [%  E     Area     Density    G     Inertia     a        b      t       
          E_s,    0,     rho_s,    G_s,       0,     a_s,    b_s,    t_s;  % Spars material  = 1
          E_r,    0,     rho_r,    G_r,       0,     a_r,    b_r,    t_r   % Ribs material   = 2
          ]';
% Material:
%   - materiall 1 = mat(:, 1)
%   - materiall 2 (E_r) = mat( 1, 2 )
      


%% Fixed nodes

Fnod = [    % Node       DOF     Displacement
               1,         1,         0;
               1,         2,         0;
               1,         3,         0;
               1,         4,         0;
               1,         5,         0;
               1,         6,         0;
               5,         1,         0;
               5,         2,         0;
               5,         3,         0;
               5,         4,         0;
               5,         5,         0;
               5,         6,         0;
]';

%% Problem dimensions

Ndim = size(x,1)*2;           % DoF for each node
Nnodes = size(x,2);           % Number of nodes
Nelements = size(Tnod,2);     % Number of elements
NnodesXelement = size(Tnod,1);% Number of nodes per element
Ndofs = Ndim*Nnodes;          % Total DoF of the system


%% Connectivity matrix - T2

T2=zeros(Nelements,Ndim*NnodesXelement);
T1=Tnod'; 
p=Ndim-1;
for i=1:Nelements
    for j=1:NnodesXelement
        for e=0:(Ndim-1)
        t2(i,(j*Ndim-p+e))=Ndim*T1(i,j)-p+e;
        end
    end
end

T2=t2';


%% Aerodynamic load computation
% Ribs elements: 1 to 27
% Front spars:   28 to 53
% Rear spars:    54 to 79


Lift = zeros(Nelements,1);
Drag = zeros(Nelements,1);


F = griddedInterpolant(yy,Lift_F);
fun_lift = @(t) F(t);
F = griddedInterpolant(yy,Drag_F);
fun_drag = @(t) F(t);

%% READ LIFT AND DRAG FROM XFOIL OR SIMMILAR



for e = 1:Nelements
    
    ye1 = x(2,Tnod(1,e));     % Node 1
    ye2 = x(2,Tnod(2,e));     % Node 2 
    y_el = (ye1 + ye2)/2;

    if e >= 28 && e <= 53
       [val1,idx1]=min(abs(yy-ye1));
        [val2,idx2]=min(abs(yy-ye2));
        drag_int = integral(fun_drag, yy(idx1), yy(idx2));
        Drag(e) = drag_int;
        lift_int = integral(fun_lift, yy(idx1), yy(idx2));
        Lift(e) = lift_int;
        
        
    end
    
    if  e>=54 && e<=79
        [val1,idx1]=min(abs(yy-ye1));
        [val2,idx2]=min(abs(yy-ye2));
        lift_int = integral(fun_lift, yy(idx1), yy(idx2));
        Lift(e) = lift_int;
    end
    
end

if a == 2
    Lift(:) = 0;
    Drag(:) = 0;
end
if a == 3
    Lift = Lift*1.2;
    Drag = Drag*1.2;
end
if a == 4
    Lift = Lift*0.8;
    Drag = Drag*1.2;
end
% 
% for e = 1:Nelements
%     
%     % Nodal position y
%     ye1 = x(2,Tnod(1,e));     % Node 1
%     ye2 = x(2,Tnod(2,e));     % Node 2 
%     
%     y_el = (ye1 + ye2)/2;
%     % Front spar
%     if e >= 28 && e <= 53
%         y_1=x(2,1); 
%         y_2=x(2,29); 
%         y_3=x(2,53); 
%         % Drag Compute
%         d_1 = 1e3;   % [N/m]
%         Drag(e) = d_1*(1 - ((y_el-y_1)/(y_3 - y_1))^2);
%         % Lift Compute
%         l_1 = 15e3;  % [N/m]
%         l_2 = 4.5e3; % [N/m]
%         if y_el > y_1 && y_el < y_2
%             % Formula 1
%             Lift(e) = 1/2* (l_1 + l_2 + (l_1-l_2)* cos(pi*((y_el-y_1)/(y_2-y_1))));
%         else
%             % Formula 2
%             Lift(e) = l_2*cos((pi/2)*((y_el-y_2)/(y_3-y_2)));
%         end  
%     end
%     
%     % Rear spar
%     if  e>=54 && e<=79
%         y_1=x(2,5); 
%         y_2=x(2,32); 
%         y_3=x(2,54); 
%         % Lift Compute 
%         l_1 = 4.5e3;  % [N/m]
%         l_2 = 4.5e3; % [N/m]
%         if y_el > y_1 && y_el < y_2
%             % Formula 1
%             Lift(e) = 1/2* (l_1 + l_2 + (l_1-l_2)* cos(pi*((y_el-y_1)/(y_2-y_1))));
%         else
%             % Formula 2
%             Lift(e) = l_2*cos((pi/2)*((y_el-y_2)/(y_3-y_2)));
%         end  
%     end
%     
% 
%     
% end

end
