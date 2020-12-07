% AA554 HW#6
% U-g Method
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

C_h = 2*zeta_h*w_h*m_hh;
C_alpha = 2*zeta_alpha*w_alpha*m_alphaalpha;
C_Beta = 2*zeta_Beta*w_Beta*m_BetaBeta;
C = [C_h,0,0;0,C_alpha,0;0,0,C_Beta];

G_str = C*inv(K); 
% g_str = 2*zeta approximation can be used

% Range of reduced frequency (1/k)
k_inverse = linspace(0.01,7,50);

% Changing hinge spring K_Beta
Percentage = 100;
K(3,3) = K(3,3)*(Percentage/100); 

% Calculate g, omega, and U
[g,omega,U] = ug_solve(a,b,c,ndof,M,K,G_str,rho,k_inverse);

g_h(:,1) = g(1,1,:);
g_alpha(:,1) = g(2,1,:);
g_Beta(:,1) = g(3,1,:);

omega_h(:,1) = (omega(1,1,:))./(2*pi);
omega_alpha(:,1) = (omega(2,1,:))./(2*pi);
omega_Beta(:,1) = (omega(3,1,:))./(2*pi);

U_h(:,1) = U(1,1,:);
U_alpha(:,1) = U(2,1,:);
U_Beta(:,1) = U(3,1,:);

% Plots
figure(1), hold all
plot(U_h,g_h), plot(U_alpha,g_alpha), plot(U_Beta,g_Beta)
axis([0,50,-6,1]), xlabel("U (m/s)"), ylabel("g")
legend("h","\alpha","\beta",'Location','Best')

figure(2), hold all
plot(U_h,omega_h), plot(U_alpha,omega_alpha), plot(U_Beta,omega_Beta)
axis([0,50,0,25]), xlabel("U (m/s)"), ylabel("\omega (Hz)"), 
legend("h","\alpha","\beta",'Location','Best')

% Find U when g=0 and omega at corresponding U
flutter_speed = U_alpha(abs(g_alpha)<0.0062)
flutter_frequency = omega_alpha(U_alpha==flutter_speed)

function [g,omega,U] = ug_solve(a,b,c,ndof,M,K,G_str,rho,k_inverse)
    T = NewTfunctions(a,c);
    K_bar = (eye(ndof) + 1i.*G_str)*K;
    B_bar = -K_bar;
    
    k = 1./k_inverse;
    for i = 1:length(k)
        [FF(i),GG(i)] = NewTheodorsenFunction(k(i));
        AA(:,:,i) = NewAeroMatrix(k(i),a,c,b,T,FF(i),GG(i),ndof);
        AA_bar(:,:,i) = -1*M-((0.5*rho*b^2)./(k(i).^2)).*AA(:,:,i);
        [V(:,:,i),D(:,:,i)] = eig(AA_bar(:,:,i)/(B_bar));
        Lambda(:,:,i) = D(:,:,i)*ones(3,1);
        g(:,:,i) = (imag(Lambda(:,:,i))./real(Lambda(:,:,i)));
        omega(:,:,i) = abs(sqrt(1./real(Lambda(:,:,i))));
        U(:,:,i) = (omega(:,:,i).*b)./k(i);
    end
end