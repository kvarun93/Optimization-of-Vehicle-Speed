

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Digital Control of Automotive Powertrains - Spring 2019
%
% Model of I.C. engine dynamics for idle speed control.
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
% clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load engine geometric parameters and constant inputs

Vd = 2.4e-3;   % Displacement (m^3)
Z = 4;         % Number of Cylinders
Vm = 5.8e-3;   % Intake Manifold Volume (m^3)
J = 0.0789;    % Mass moment of inertia

p_amb = 1.0121*1e5;
T_amb = 302;
R=288;
gam = 1.35;

P0 = 26431;   % Initial MAP
N0 = 828;     % Initial RPM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model parameters (from steady-state calibration)

a = [1.69e-07,-1.136e-06,6.89e-06];  % Throttle
si = 0.812;   yi = 0.0633;           % Volumetric efficiency
P_IMEP = [0.0220,-1.1649];           % IMEP
par_sp = [-0.0017 -0.0277 1.36];     % Spark timing effect
par_fr = [7.4198e-7 -4.989e-4 11.3]; % Friction
par_IMEP0 = [1.2323e-4 2.1256];      % Base IMEP model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversion to Crank Angle Domain

% Coefficients for crank-angle based model
Kthr = p_amb/sqrt(R*T_amb)*sqrt(gam)*sqrt((2/(gam+1))^((gam+1)/(gam-1)));
Kp1 = R*T_amb/Vm;
Kp2 = si*Vd/(4*pi*Vm);
Kp3 = yi*Vd/(4*pi*Vm);
KT = 1e5*Vd/(4*pi);
Kfr1 = (30/pi)^2 * par_fr(1);
Kfr2 = (30/pi) * par_fr(2);
Kfr3 = par_fr(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Equilibrium Condition (p_0 T_0 w_0)

setp(1) = 9.81; % Throttle Position
setp(2) = -25;  % Spark Timing
setp(3) = 10;   % Load Torque
X0 = [26424 21.3765773202354 83.9019428270409]; % Equilibrium Conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearization

% Coefficients for linearized model as shown in lecture
K1 = Kp1*Kthr*(2*a(1)*setp(1)+a(2));
K2 = Kp1*Kthr*(a(1)*setp(1)^2 + a(2)*setp(1) + a(3));
K3 = Kp2;
Kp = KT*par_IMEP0(1)*(par_sp(1)*setp(2)^2 + par_sp(2)*setp(2) + par_sp(3));    % Pressure multiplier
Kt = KT*(par_IMEP0(1)*X0(1) - par_IMEP0(2)) * (par_sp(1)*setp(2) + par_sp(2)); % Spark Timing multiplier
Kf = 2*Kfr1*X0(3)^2 + Kfr2*X0(3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants Calculation
Tau_d = pi/3;
Ndel = 3; 

CT1 = 1+ ((pi*K3)/6);
CT2 = ((pi*K3)/6)-1;
CT3 = (pi*K2)/(3*2*X0(3)^2);
CT4 = (pi*K1)/(2*X0(3)*3);
CT5 = ((pi*Kf)/(2*3*J*X0(3)^2))+1;
CT6 = ((pi*Kf)/(2*J*(X0(3)^2)*3))-1;
CT7 = pi/(3*2*J*X0(3));
CT8 = pi/(3*2*J*X0(3));


%%

% ------ Closed Loop (Hand Tuning) -------- %

KP = 0.8;   
TI = 1; 
KI = 0.005;
Td = 1;   
Kd = 4


%%

% ------ Closed Loop (fmincon-contrained optimization) -------- %

x0 = [0.5,0.005,4]; % [0.01,0.001,1] % Setting the inital gains
lb = [0 0 0]; % the lower bounds 
ub = [2 2 2]; % Setting upper bounds
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp'); % fmincon using sqp algorithm 
optionsT = optimoptions(@fmincon,'Display','iter','FunctionTolerance',1e-10) % changing the Tolerance value
options = optimoptions('fmincon','Display','iter','Algorithm','active-set'); % Trying out with different algorithm
optionsP = optimoptions('patternsearch','Display','iter'); % Pattern search algorithm
[x,fval,exitflag,output] = fmincon(@myFunction,x0,A,b,Aeq,beq,lb,ub,nonlcon,options); % Problem 2 (fmincon)
x = fmincon(@myFunction,x0,A,b,Aeq,beq,lb,ub,nonlcon,optionsT); % Problem 3 fmincon for different tolerance level 
x = patternsearch(@myFunction,x0,A,b,Aeq,beq,lb,ub,nonlcon,optionsP) % Problem 3 pattern search is suggested

%%

% ------ Closed Loop (fminsearch-uncontrained optimization) -------- %

 x0 = [1,1,1]; % Setting the initial gains
 options = optimset('Display','iter') 
 [x,fval,exitflag,output] = fminsearch(@myFunction,x0,options) 

%%

% ------ Closed Loop (fminunc-uncontrained optimization) -------- %

x0 = [0.01,0.001,2]; 
options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton'); % With quasi-Newton Algorithm
[x,fval,exitflag,output] = fminunc(@myFunction,x0,options)


%%

% ------ General Function Definition -------- %
% 
function F = myFunction(x)
    assignin('base','KP',x(1))
    assignin('base','KI',x(2))
    assignin('base','Kd',x(3))
    simout = sim('idle_speed_model_project_pid');
    
  % ------ summation of absolute value speed error objective Function -------- %
  
    E = integral(@(t)(abs(error.signals.values)),0,numel(error.time),'ArrayValued',true);
    F = sum(E) % This is function for summation for absolute value speed error.

  % ------ Norm-squared speed objective Function -------- % 
  
    E = integral(@(t)power((norm(error.signals.values)),2),0,numel(error.time),'ArrayValued',true)
    F = sum(E)
  
% ------ Time((error)^2) ocjective Function  -------- %

    E = integral(@(t)power(error.signals.values,2).*(error.time),0,numel(error.time),'ArrayValued',true)
    F = sum(E)   
   
end
 