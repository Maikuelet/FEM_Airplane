function [yy,Lift,Drag] = WingData_function() 

%% INPUT DATA%% INPUT DATA

%Wind data
v = 221; % Reference velocity [m/s]
h = 1; % Height     [m]
Surf_FS  = 122.5; % Flight stream reference surface [ m^2]
v_FS = 221;        % FS ref velocity [m/s]
L_FS = 624932; % FS Lift
MTOW = 70000;
% ISA ref param
R   = 287;              % Ideal gas constant [J kg^-1 K^-1]
Gamma = 1.4;            % Diatomic gas ct (Air)

a = -6.6/1000;          % Height ct [ºC/m]
T_ref = 273.15;         % Temp ref [K]
T_0 = T_ref + 15;       % Temp SL  [K]
rho_0 = 1.225;          % Density SL [kg/m^3]
g_0 = 9.81;

T = T_0 + a*h;
rho = rho_0*exp((-g_0*h)/(R*T)); 

rho = 0.366;

%% DATA LOADING
%Define File names and Airfoil Name
pathh     = pwd;
myfolder  = 'Data';

Name1 = 'cd_wing_cruise.txt';
Name2 = 'cl_distrib_read3.txt';
Name3 = 'WingCoord.xlsx';


% Drag Distribution
f1 = fullfile(pathh , myfolder, Name1);
D = dlmread(f1);

% Lift Distribution
f2 = fullfile(pathh , myfolder, Name2);
% fileID = fopen(f2);
% L = fscanf(fileID,%f);
L = dlmread(f2);

% Wing Coordinates
f3 = fullfile(pathh , myfolder, Name3);
WC_LE = xlsread(f3,'LE')/1000;
WC_TE = xlsread(f3,'TE')/1000;




%% PREPROCESS

% LE TE Spline & Pchip generation
P = 0.001;  %Spline Parameter
mode = 0;       % Mode for spline or Pchip

yy = 0:P:WC_LE(size(WC_LE,1),2);

if mode == 1
S_LE = spline(WC_LE(:,2),WC_LE(:,1),yy);
S_TE = spline(WC_TE(:,2),WC_TE(:,1),yy);
else
S_LE = pchip(WC_LE(:,2),WC_LE(:,1),yy);
S_TE = pchip(WC_TE(:,2),WC_TE(:,1),yy);
end

% Chord calculation

N = size(S_LE,2);

chord = zeros(1,N);
chord(:) = S_TE(:) - S_LE(:);
S = sum(chord(:));

% Drag treatment
k = size(D,1);
cd = D(:,2);
y_cd = D(:,1);
diff = y_cd(k) - WC_LE(size(WC_LE,1),2);
y_cd(2:end) = y_cd(2:end) - diff;

S_cd = pchip(y_cd,cd,yy);


F = griddedInterpolant(yy,S_cd);
fun_cd = @(t) F(t);
cd_int = integral(fun_cd, yy(1), yy(end));

Cd_2Wings = cd_int*2;


Drag = (S_cd.*chord*rho*v^2)/2;     % One wing drag

% Integral Lift
F = griddedInterpolant(yy,Drag);
fun_drag = @(t) F(t);
drag_int = integral(fun_drag, yy(1), yy(end));

% Lift Treatment
k = size(L,1);
i=1;
j=1;


for ii=1:k
    if L(ii,6) > 0
        L_Up(i) = L(ii,6);
        y_Up(i) = L(ii,5);
        i = i+1;
    end
    if  L(ii,6) < 0
        L_Dn(j) = L(ii,6);
        y_Dn(j) = L(ii,5);
        j= j +1;
    end
end

L_Up(1) = [];
y_Up(1) = [];
kU = size(L_Up,2);
L_Up(kU) = [];
y_Up(kU) = [];
kU = size(L_Up,2);
L_Up(kU/2+1) = [];
y_Up(kU/2+1) = [];
kU = size(L_Up,2);
L_Up(floor(kU/2)+1) = [];
y_Up(floor(kU/2)+1) = [];

kD = size(L_Dn,2);
L_Dn(kD/2+1) = [];
y_Dn(kD/2+1) = [];
kD = size(L_Dn,2);
L_Dn(floor(kD/2)+1) = [];
y_Dn(floor(kD/2)+1) = [];

kU = size(L_Up,2);
b_wing = abs(y_Up(1) - y_Up(kU/2));  % Winglegth from Flighstream

    
cl_Fs = L_Up +L_Dn;     % C_L sum up and down

k = size(cl_Fs,2);
j=1;
for i = 1:k
   if  y_Up(i) > 0
       cl_vec(j) = cl_Fs(i);
       y_vec(j)  = y_Up(i);
       j = j+1;
    end
end

y_vec = y_vec - y_vec(1);
k = size(y_vec,2);
y_vec(k-1) = y_vec(k-1) + (WC_LE(size(WC_LE,1),2) - y_vec(k))/2;
y_vec(k) = WC_LE(size(WC_LE,1),2); % Change last postion to match Catia model

S_cl = pchip(y_vec,cl_vec,yy);  %spline(y_vec,cl_vec,yy);


% Integral S_cl

F = griddedInterpolant(yy,S_cl);
fun_cl = @(t) F(t);
cl_int = integral(fun_cl, yy(1), yy(end));

Cl_2Wings = cl_int*2;

% Lift Computation !!!!!!!!!!!!!!!!!!!!!!!


Lift = (S_cl.*chord*rho*v^2)/2;     % One wing lift

% Integral Lift
F = griddedInterpolant(yy,Lift);
fun_lift = @(t) F(t);
lift_int = integral(fun_lift, yy(1), yy(end));

lift_int2 = 0;
for i =1:(size(yy,2)-1)
    
    lift_int2 = lift_int2 + Lift(i)*(yy(i+1)-yy(i));

end

% Integral Wing surface
F = griddedInterpolant(yy,chord);
fun_Surf = @(t) F(t);
S_int = integral(fun_Surf, yy(1), yy(end));
 
x = (MTOW*g_0)/lift_int;
Lift = Lift*x;
Drag = Drag*3;
end