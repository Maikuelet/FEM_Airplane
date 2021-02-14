function  PostProcess(M_struct,M_rib,M_spar,V_rib,V_spar,Lift_T,Drag_T,V_tot,U,u_int,F_int,El_L,N_x,Q_y,Q_z,T_x,M_y,M_z,R,Rr,El_W)

InputData;

% Plot tittles with lattex typography
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultaxesticklabelinterpreter','latex');
% set(groot,'defaultlegendinterpreter','latex');


%% DATA DISPLAY - PLOTS

% PLOT 1 - PLOTWING
x = x.';
Tnod = Tnod.';
U = U.';
plotWing(x,Tnod,El_L,U,u_int,N_x,Q_y,Q_z,T_x,M_y,M_z);

% PLOT 2 - PLOTBEAMS3D DEFORMATION
nsub  = 10;     % Number of subdivisions
factor = 2.5;   % Amplification factor

plotBeams3D_def(x,Tnod,nsub,El_L,u_int,factor,R)

% PLOT BEAMS 3D - CRITICAL BEAMS
maxAxial_e = find(N_x(1,:) == max(N_x));

maxShear_y = max(abs(Q_y));
maxShear_z = max(abs(Q_z));

if maxShear_y > maxShear_z
    maxShear_e = find(abs(Q_y(1,:)) == maxShear_y);
    maxShear = maxShear_y;
    dirs = 'Y';
else
    maxShear_e = find(abs(Q_z(1,:)) == maxShear_z);
    maxShear = maxShear_z;
    dirs = 'Z';
end

maxBend_y = max(max(abs(M_y)));
maxBend_z = max(max(abs(M_z)));

if maxBend_y > maxBend_z
    maxBend_e = find(abs(M_y(1,:)) == maxBend_y);
    maxBend = maxBend_y;
    dirb = 'Y';
else
    maxBend_e = find(abs(M_z(1,:)) == maxBend_z);
    maxBend = maxBend_z;
    dirb = 'Z';
end


e_plot = maxAxial_e;
plotBeams3D(Tnod,nsub,El_L,u_int,e_plot)


%% DATA DISPLAY - CONSOLE
j=0;
for i =1:12
    if i<=6
        fprintf('Node:  1  Dof: %d  Reaction: %0.4f  \n',i,Rr(i));
    else
        j=j+1;
        fprintf('Node:  5  Dof: %d  Reaction: %0.4f  \n',j,Rr(i));    
    end    
end
disp('----------------------------------------------------------------')
fprintf('Structure Mass %0.4f  Kg \n',M_struct);
fprintf('Structure Volume %0.4f  m^3 \n',V_tot);

fprintf('Ribs Mass %0.4f  Kg \n',M_rib);
fprintf('Ribs Volume %0.4f  m^3 \n',V_rib);

fprintf('Spars Mass %0.4f  Kg \n',M_spar);
fprintf('Spars Volume %0.4f  m^3 \n',V_spar);

fprintf('Wing - Struct Mass %0.4f  Kg \n',-M_struct+M_w);
fprintf('Element Mass %0.4f  Kg \n',sum(El_W)/g);
fprintf('Wing & Motor Mass %0.4f  Kg \n',M_struct+M_w+M_eng);
fprintf('Wing Lift %0.4d  N \n',Lift_T);
fprintf('Wing Drag %0.4f  N \n',Drag_T);
disp('----------------------------------------------------------------')

fprintf('Max Axial force = %d    At element:  %d   \n',max(N_x),maxAxial_e);
fprintf('Max Shear force = %d    At element:  %d  Direction: %s \n',maxShear,maxShear_e,dirs);
fprintf('Max Bending moment = %d    At element:  %d  Direction: %s  \n',maxBend,maxBend_e,dirb);
disp('----------------------------------------------------------------')


end