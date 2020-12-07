% 3DOF 2D Wing+Control Surface 
clear all; close all; clc

% Inputs
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


% Mass matrix (using actual masses then dividing by span)
M = [(m+m_blocks)/span, S_alpha, S_Beta;
      S_alpha, I_alpha, I_Beta+S_Beta*b*(c-a);
      S_Beta, I_Beta+S_Beta*b*(c-a), I_Beta]
% Stiffness matrix
K = [K_h, 0, 0;
     0, K_alpha, 0;
     0, 0, K_Beta]
 
 
% K*Phi = Lambda*M*Phi
[Phi,Lambda] = eig(K,M);
Omega = sqrt(Lambda); % natural frequencies in rad/s
Frequency = Omega./(2*pi) % natural frequencies in Hz
% Measured values of the natural frequencies of the coupled system:
% 4.375, 9.125, and 18.625 Hz (27.489, 57.334, 117.024 rad/s)


% Plot
% Take the origin to be at the elastic axis
h_max = max(Phi(1,:));
h_min = min(Phi(1,:));
for i = 1:3
    Phi_norm = Phi(:,i)/norm(Phi(:,i));
    h = Phi_norm(1);
    alpha = Phi_norm(2);
    Beta = Phi_norm(3);
    P1x = -(b+b*a)*cos(alpha);
    P1y = (b+b*a)*sin(alpha)-h;
    P2x = b*(c-a)*cos(alpha);
    P2y = -b*(c-a)*sin(alpha)-h;
    P3x = P2x+(b-b*c)*cos(alpha+Beta);
    P3y = P2y-(b-b*c)*sin(alpha+Beta);
    figure()
    plot([-(b+b*a), b-(b*a)], [0, 0],'k-','LineWidth',2), hold on
    plot([P1x,P2x], [P1y,P2y],'b-','LineWidth',2), hold on
    plot([P2x,P3x], [P2y,P3y],'r-','LineWidth',2)
    title(['Mode Shape ', num2str(i)]), xlabel('x(m)'), ylabel('y(m)') 
    axis([-2*b 2*b, -(2*b)-h_max 2*b-h_min]), legend({'baseline','airfoil','aileron'},'Location','Best')
    daspect([1 1 1])
    %d_wing = sqrt((P2x-P1x)^2+(P2y-P1y)^2)
    %d_aileron = sqrt((P3x-P2x)^2+(P3y-P2y)^2)
end
