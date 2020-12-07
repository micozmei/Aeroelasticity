% AA554 HW#8
% Control
clear all; close all; clc

% Inputs
ndof = 3;             %number of degrees of freedom
rho = 1.225;          %density of air
chord = 0.254;        %chord length in meters
b = chord/2;          %mid-chord length in meters
a = -0.5;             %elastic axis location behind mid-chord (normalized w.r.t. b)
c = 0.5;              %hinge line location behind mid-chord (normalized w.r.t. b)
span = 0.52;          %wing span in meters
m_wing = 0.62868;     %mass of wing only (no aileron) in kg (actual mass before dividing by span)
m_Beta = 0.18597;     %mass of aileron in kg (actual mass before dividing by span)
m = m_wing + m_Beta;  %total mass of wing and aileron in kg (actual mass before dividing by span)
m_blocks = 2*0.47485; %mass of support blocks (affecting only plunge motion) in kg (before dividing by span)
m_span = 1.558;       %mass/span of wing-aileron in kg/m (note, this is not consistent with the wing and aileron masses)
S_alpha = 0.08587;    %static unbalance of whole airfoil-aileron about its elastic axis in kg/m (per unit span)
S_Beta = 0.00395;     %static unbalance of aileron on its hinge in kg/m (per unit span)
I_alpha = 0.01347;    %pitch moment of inertia of configuration about its elastic axis in kg^2/m (per unit span)
I_Beta = 0.0003264;   %pitch moment of inertia of aileron about its hinge in kg^2/m (per unit span)
K_h = 2818.4;         %spring stiffness in Newton/m per unit span (Newton/m^2)
K_alpha = 37.34;      %spring stiffness Newton*m/rad per unit span (Newton/rad)
K_Beta = 3.895;       %spring stiffness Newton*m/rad per unit span (Newton/rad)
zeta_h = 0.0113;      %given damping ratios
zeta_alpha = 0.01626;
zeta_Beta = 0.0115;

% Mass matrix (using actual masses then dividing by span)
M = [(m+m_blocks)/span, S_alpha, S_Beta;
     S_alpha, I_alpha, I_Beta+S_Beta*b*(c-a);
     S_Beta, I_Beta+S_Beta*b*(c-a), I_Beta];

m_hh = M(1,1);
m_alphaalpha = M(2,2);
m_BetaBeta = M(3,3);
  
% Stiffness matrix
K = [K_h, 0, 0;
     0, K_alpha, 0;
     0, 0, K_Beta];

% Structural damping matrix
w_h = sqrt(K_h/m_hh);
w_alpha = sqrt(K_alpha/m_alphaalpha);
w_Beta = sqrt(K_Beta/m_BetaBeta);
w = [w_h,w_alpha,w_Beta];

C_h = 2*zeta_h*w_h*m_hh;
C_alpha = 2*zeta_alpha*w_alpha*m_alphaalpha;
C_Beta = 2*zeta_Beta*w_Beta*m_BetaBeta;
C = [C_h,0,0;0,C_alpha,0;0,0,C_Beta];

% Aerodynamic matrix
k = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9,1.1,1.4,1.7,2.0,2.5,3.0];
AA = zeros(ndof,ndof,length(k));
for i = 1:length(k)
   [FF,GG] = NewTheodorsenFunction(k(i));
   [T] = NewTfunctions(a,c);
   [AA(:,:,i)] = NewAeroMatrix(k(i),a,c,b,T,FF,GG,ndof);
end

% Aerodynamic lag terms
Beta_bar_1 = 0.3;
Beta_bar_2 = 0.8;

% Left hand side matrix
Left = zeros(2*(length(k)-1),4);
for i = 1:length(k)-1
    Left(2*i-1,:) = [0, -k(i+1)^2, k(i+1)^2/(k(i+1)^2+Beta_bar_1^2), k(i+1)^2/(k(i+1)^2+Beta_bar_2^2)];
    Left(2*i,:) = [k(i+1), 0, k(i+1)*Beta_bar_1/(k(i+1)^2+Beta_bar_1^2), k(i+1)*Beta_bar_2/(k(i+1)^2+Beta_bar_2^2)];
end

% Right hand side matrix
Right = zeros(2*(length(k)-1),1,ndof^2);
count = 0;
for m = 1:ndof
    for n = 1:ndof
        count = count+1;
        for i = 1:length(k)-1
            Right(2*i-1,1,count) = real(AA(m,n,i+1))-real(AA(m,n,1));
            Right(2*i,1,count) = imag(AA(m,n,i+1));
        end
    end
end

% P_bar matrix
P_bar = zeros(4,1,ndof^2);
for i = 1:ndof^2
    P_bar(:,:,i) = Left\Right(:,:,i);
end

% Roger Approximation with fine mesh of k values
k_fine = linspace(0,3,101);
AA_Roger = zeros(ndof,ndof,length(k_fine));
for i = 1:length(k_fine)
    count = 0;
    for m = 1:ndof
        for n = 1:ndof
            count = count+1;
            AA_Roger(m,n,i) = AA(m,n,1)+1i*k_fine(i)*P_bar(1,1,count)-k_fine(i)^2*P_bar(2,1,count)+...
                (k_fine(i)^2+1i*k_fine(i)*Beta_bar_1)/(k_fine(i)^2+Beta_bar_1^2)*P_bar(3,1,count)+...
                (k_fine(i)^2+1i*k_fine(i)*Beta_bar_2)/(k_fine(i)^2+Beta_bar_2^2)*P_bar(4,1,count);
        end
    end
end
 
% State Space & Root Locus
k_max = 2;
w_max = 2*round(max(w));
U_min = w_max*b/k_max;
U_flutter = 22.1208; %estimated flutter speed
U = U_flutter*1.2; %linspace(U_min,30,1000);
% test factors of 0.8, 0.9, 1.0, 1.1, and 1.3 as well

% P_bar matrix components
P0_bar = AA(:,:,1);
P1_bar = zeros(ndof);
P2_bar = zeros(ndof);
P3_bar = zeros(ndof);
P4_bar = zeros(ndof);
count = 0;
for m = 1:ndof
    for n = 1:ndof
        count = count+1;
        P1_bar(m,n) = P_bar(1,1,count);
        P2_bar(m,n) = P_bar(2,1,count);
        P3_bar(m,n) = P_bar(3,1,count);
        P4_bar(m,n) = P_bar(4,1,count);
    end
end

% Proportional Gains (positive feedback)
P_h = 0;
P_alpha = 0;
P_matrix = [0,0,0; 0,0,0; P_h,P_alpha,0];
R_h = 0;
R_alpha = 0;
R_matrix = [0,0,0; 0,0,0; R_h,R_alpha,0];

% Calculate Eigenvlaues
gain_sweep = 0; %use this to test multiple gain values -> linspace(-1,1,100)
Lambda = zeros(4*ndof,length(gain_sweep));

for i = 1:length(gain_sweep)
    % R_matrix = [0,0,0; 0,0,0; R_h,gain_sweep(i),0];%use this to test multiple gain values
    P0 = P0_bar;
    P1 = P1_bar*(b/U);
    P2 = P2_bar*(b/U)^2;
    P3 = P3_bar;
    P4 = P4_bar;
    Beta_1 = Beta_bar_1*(U/b);
    Beta_2 = Beta_bar_2*(U/b);
    M_bar = M-(0.5*rho*U^2)*P2; % M_bar = M_bar_bar
    C_bar = C-(0.5*rho*U^2)*P1;
    C_bar_bar = C_bar-R_matrix;
    K_bar = K-(0.5*rho*U^2)*P0;
    K_bar_bar = K_bar-P_matrix;
    s = [zeros(ndof), eye(ndof), zeros(ndof), zeros(ndof);
         -inv(M_bar)*K_bar_bar, -inv(M_bar)*C_bar_bar, (0.5*rho*U^2)*inv(M_bar)*P3, (0.5*rho*U^2)*inv(M_bar)*P4;
         zeros(ndof), eye(ndof), -Beta_1*eye(ndof), zeros(ndof);
         zeros(ndof), eye(ndof), zeros(ndof), -Beta_2*eye(ndof)];
    Lambda(:,i) = reshape(eig(s),[],1);
end

% Root Locus Plot
figure(1) 
max_root = -inf;
for i = 1:4*ndof
    plot(real(Lambda(i,:)),imag(Lambda(i,:)),'*'), hold on
    if max_root < real(Lambda(i,:))
        max_root = real(Lambda(i,:));
    end
end
xlabel('\sigma'), ylabel('j \omega'), title('Root Locus')
legend('location', 'NorthWest')
% max_root

% Linear Quadratic Regulator
A_matrix = s;
B_matrix = [0;0;0;
            inv(M_bar)*[0;0;1];
            0;0;0;
            0;0;0];
C_matrix = [1,0,0,0,0,0,0,0,0,0,0,0;
            0,1,0,0,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0];
D_matrix = 0;
sys = ss(A_matrix,B_matrix,C_matrix,D_matrix);

Q = C_matrix'*C_matrix;
R = 1; % Penalize actuator input
[G,S,e] = lqr(A_matrix,B_matrix,Q,R); % G = optimal gain matrix where {u}=-[G]{x}
sys_cl = ss(A_matrix-B_matrix*G,B_matrix,C_matrix,D_matrix);

figure(2)
for i = 1:length(e) % e = eig(A-B*G)
    plot(real(e(i)),imag(e(i)),'*'), hold on
end
xlabel('\sigma'), ylabel('j \omega'), title('Closed Loop Poles')
legend('location', 'NorthWest')

% Plot Time History for sys (open loop) and sys_cl (closed loop)
t = 0:0.001:5;
u = zeros(length(t),1);
x0 = [1,0,0,0,0,0,0,0,0,0,0,0]; % Unit plunge at t=0
figure(3)
lsim(sys,u,t,x0), title('Open Loop Time History')
figure(4)
lsim(sys_cl,u,t,x0), title('Closed Loop Time History')


%% Gain tuned for 1.2*U_flutter used for other speeds
U_test = [0.8,0.9,1.0,1.1,1.2,1.3];
U = U_test*U_flutter;

poles = zeros(ndof*4,length(U));
for u = 1:length(U)
    Lambda = zeros(4*ndof,length(gain_sweep));
    for i = 1:length(gain_sweep)
        % R_matrix = [0,0,0; 0,0,0; R_h,gain_sweep(i),0];%use this to test multiple gain values
        P0 = P0_bar;
        P1 = P1_bar*(b/U(u));
        P2 = P2_bar*(b/U(u))^2;
        P3 = P3_bar;
        P4 = P4_bar;
        Beta_1 = Beta_bar_1*(U(u)/b);
        Beta_2 = Beta_bar_2*(U(u)/b);
        M_bar = M-(0.5*rho*U(u)^2)*P2; % M_bar = M_bar_bar
        C_bar = C-(0.5*rho*U(u)^2)*P1;
        C_bar_bar = C_bar-R_matrix;
        K_bar = K-(0.5*rho*U(u)^2)*P0;
        K_bar_bar = K_bar-P_matrix;
        s = [zeros(ndof), eye(ndof), zeros(ndof), zeros(ndof);
             -inv(M_bar)*K_bar_bar, -inv(M_bar)*C_bar_bar, (0.5*rho*U(u)^2)*inv(M_bar)*P3, (0.5*rho*U(u)^2)*inv(M_bar)*P4;
             zeros(ndof), eye(ndof), -Beta_1*eye(ndof), zeros(ndof);
             zeros(ndof), eye(ndof), zeros(ndof), -Beta_2*eye(ndof)];
        Lambda(:,i) = reshape(eig(s),[],1);
    end
    A_matrix = s;
    B_matrix = [0;0;0;
                inv(M_bar)*[0;0;1];
                0;0;0;
                0;0;0];
    sys = ss(A_matrix,B_matrix,C_matrix,D_matrix);
    sys_cl = ss(A_matrix-B_matrix*G,B_matrix,C_matrix,D_matrix);
    poles(:,u) = pole(sys_cl);
end

figure(5)
for i = 1:length(U)
    plot(poles(:,i),'*', 'displayname', strcat('U_{flutter}*',num2str(U_test(i)))), hold on
end
xlabel('\sigma'), ylabel('j \omega'), title('Closed Loop Poles')
legend('location', 'NorthWest')