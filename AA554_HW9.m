% AA554 HW#9
% Gusts on 2D & 2 DOF Wing
clear all; close all; clc

% Inputs
ndof = 2;             %number of degrees of freedom
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
zeta_alpha = 0.01626; %given damping ratios
zeta_Beta = 0.0115;   %given damping ratios

% Mass matrix (using actual masses then dividing by span)
M = [(m+m_blocks)/span, S_alpha;
     S_alpha, I_alpha];

m_hh = M(1,1);
m_alphaalpha = M(2,2);
  
% Stiffness matrix
K = [K_h, 0;
     0, K_alpha];

% Structural damping matrix
w_h = sqrt(K_h/m_hh);
w_alpha = sqrt(K_alpha/m_alphaalpha);
w = [w_h,w_alpha];

C_h = 2*zeta_h*w_h*m_hh;
C_alpha = 2*zeta_alpha*w_alpha*m_alphaalpha;
C = [C_h,0;
     0,C_alpha];
G_str = C*inv(K); % G_str = 2*zeta approximation can be used

% Aerodynamic matrix
k_aero = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9,1.1,1.4,1.7,2.0,2.5,3.0];
AA = zeros(ndof,ndof,length(k_aero));
for i = 1:length(k_aero)
   [FF,GG] = NewTheodorsenFunction(k_aero(i));
   [T] = NewTfunctions(a,c);
   [AA(:,:,i)] = NewAeroMatrix2DOF(k_aero(i),a,c,b,T,FF,GG,ndof);
end

%% Rational Function Approximation of Aerodynamic Matrix
% Aerodynamic lag terms
P_size = 4; % depends on max # of aero lag terms you want to use
Beta_bar_1 = 0.3;
Beta_bar_2 = 0.8;

% Left hand side matrix
Left = zeros(2*(length(k_aero)-1),P_size); % used P_size-1 since P5_bar=0
for i = 1:length(k_aero)-1
    Left(2*i-1,:) = [0, -k_aero(i+1)^2, k_aero(i+1)^2/(k_aero(i+1)^2+Beta_bar_1^2), k_aero(i+1)^2/(k_aero(i+1)^2+Beta_bar_2^2)];
    Left(2*i,:) = [k_aero(i+1), 0, k_aero(i+1)*Beta_bar_1/(k_aero(i+1)^2+Beta_bar_1^2), k_aero(i+1)*Beta_bar_2/(k_aero(i+1)^2+Beta_bar_2^2)];
end

% Right hand side matrix
Right = zeros(2*(length(k_aero)-1),1,ndof^2);
count = 0;
for m = 1:ndof
    for n = 1:ndof
        count = count+1;
        for i = 1:length(k_aero)-1
            Right(2*i-1,1,count) = real(AA(m,n,i+1))-real(AA(m,n,1));
            Right(2*i,1,count) = imag(AA(m,n,i+1));
        end
    end
end

% P_bar matrix
P_bar = zeros(P_size,1,ndof^2);
for i = 1:ndof^2
    P_bar(:,:,i) = Left\Right(:,:,i);
end
P_bar(5,1) = 0;
P0_bar = AA(:,:,1);
P1_bar = zeros(ndof);
P2_bar = zeros(ndof);
P3_bar = zeros(ndof);
P4_bar = zeros(ndof);
P5_bar = zeros(ndof);
count = 0;
for m = 1:ndof
    for n = 1:ndof
        count = count+1;
        P1_bar(m,n) = P_bar(1,1,count);
        P2_bar(m,n) = P_bar(2,1,count);
        P3_bar(m,n) = P_bar(3,1,count);
        P4_bar(m,n) = P_bar(4,1,count);
        P5_bar(m,n) = P_bar(5,1,count);
    end
end

% Roger Approximation with fine mesh of k values
k_fine = linspace(0,3,301);
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

%% Gust Force Vector
k = linspace(0,3,301); % increase k max value to see complete spiral solution
[F_gust, F_gust_modified, H, Phi] = Gust_vector(k,a,b);

% Plot Sears/Modified Sears Functions and Gust Force Vector
figure()
plot(real(Phi),imag(Phi),real(H),imag(H))
xlabel('Real'),ylabel('Imaginary'),title('Sears Function')
legend('Sears Function (\phi)','Modified Sears Function (H)')

figure()
plot(real(F_gust(1,:)),imag(F_gust(1,:)),real(F_gust(2,:)),imag(F_gust(2,:)))
xlabel('Real'),ylabel('Imaginary'),title('Gust Force Vector')
legend('h component','\alpha component')

figure()
plot(real(F_gust_modified(1,:)),imag(F_gust_modified(1,:)),real(F_gust_modified(2,:)),imag(F_gust_modified(2,:)))
xlabel('Real'),ylabel('Imaginary'),title('Modified Gust Force Vector')
legend('h component','\alpha component')

%% Rational Function Approximation of Gust Force Vector 
[F_gust_0, F_gust_modified_0, H_0, Phi_0] = Gust_vector(k_aero,a,b);

% Aerodynamic lag terms
P_size_gust = P_size; % depends on how many aero lag terms you want to use
Beta_bar_1_gust = 0.3;
Beta_bar_2_gust = 0.35;
Beta_bar_3_gust = 0.8;

% Left hand side matrix
Left_gust = zeros(2*(length(k_aero)-1),P_size_gust);
for i = 1:length(k_aero)-1
    Left_gust(2*i-1,:) = [0, k_aero(i+1)^2/(k_aero(i+1)^2+Beta_bar_1_gust^2), ...
        k_aero(i+1)^2/(k_aero(i+1)^2+Beta_bar_2_gust^2), k_aero(i+1)^2/(k_aero(i+1)^2+Beta_bar_3_gust^2)];
    Left_gust(2*i,:) = [k_aero(i+1), k_aero(i+1)*Beta_bar_1_gust/(k_aero(i+1)^2+Beta_bar_1_gust^2), ...
        k_aero(i+1)*Beta_bar_2_gust/(k_aero(i+1)^2+Beta_bar_2_gust^2), k_aero(i+1)*Beta_bar_3_gust/(k_aero(i+1)^2+Beta_bar_3_gust^2)];
end

% Right hand side matrix
Right_gust = zeros(2*(length(k_aero)-1),ndof);
count = 0;
for row = 1:ndof
    count = count+1;
    for i = 1:length(k_aero)-1
        Right_gust(2*i-1,count) = real(F_gust_modified_0(row,i+1))-real(F_gust_modified_0(row,1));
        Right_gust(2*i,count) = imag(F_gust_modified_0(row,i+1));
    end
end

% P_bar_gust matrix
P_bar_gust = zeros(P_size_gust,ndof);
for i = 1:ndof
    P_bar_gust(:,i) = Left_gust\Right_gust(:,i);
end
P0_bar_gust = F_gust_modified_0(:,1);
P1_bar_gust = zeros(ndof,1);
P3_bar_gust = zeros(ndof,1);
P4_bar_gust = zeros(ndof,1);
P5_bar_gust = zeros(ndof,1);
count = 0;
for i = 1:ndof
    count = count+1;
    P1_bar_gust(i) = P_bar_gust(1,count);
    P3_bar_gust(i) = P_bar_gust(2,count);
    P4_bar_gust(i) = P_bar_gust(3,count);
    P5_bar_gust(i) = P_bar_gust(4,count);
end

% Gust Force Vector
F_gust_approx = zeros(ndof,length(k_fine));
for i = 1:length(k_fine)
    count = 0;
    for row = 1:ndof
        count = count+1; % note P2 aero inertia terms not used for gusts
        F_gust_approx(row,i) = F_gust_modified_0(row,1)+1i*k_fine(i)*P_bar_gust(1,count)+...
            (k_fine(i)^2+1i*k_fine(i)*Beta_bar_1_gust)/(k_fine(i)^2+Beta_bar_1_gust^2)*P_bar_gust(2,count)+...
            (k_fine(i)^2+1i*k_fine(i)*Beta_bar_2_gust)/(k_fine(i)^2+Beta_bar_2_gust^2)*P_bar_gust(3,count)+...
            (k_fine(i)^2+1i*k_fine(i)*Beta_bar_3_gust)/(k_fine(i)^2+Beta_bar_3_gust^2)*P_bar_gust(4,count);
    end
end

figure()
plot(real(F_gust_approx(1,:)),imag(F_gust_approx(1,:)),real(F_gust_approx(2,:)),imag(F_gust_approx(2,:)))
xlabel('Real'),ylabel('Imaginary'),title('Gust Force Vector (Approximate)')
legend('h component','\alpha component')

%% System
w_max = 2*round(max(w));
U_min = w_max*b/max(k);
U = linspace(U_min,30,1000);

% Calculate eigenvlaues
Lambda = zeros(11,length(U));
for i = 1:length(U)
    
    qD = 0.5*rho*U(i)^2;
    
    P0 = P0_bar;
    P1 = P1_bar*(b/U(i));
    P2 = P2_bar*(b/U(i))^2;
    P3 = P3_bar;
    P4 = P4_bar;
    Beta_1 = Beta_bar_1*(U(i)/b);
    Beta_2 = Beta_bar_2*(U(i)/b);
    
    P0_gust = P0_bar_gust;
    P1_gust = P1_bar_gust*(b/U(i));
    P3_gust = P3_bar_gust;
    P4_gust = P4_bar_gust;
    P5_gust = P5_bar_gust;
    Beta_1_gust = Beta_bar_1_gust*(U(i)/b);
    Beta_2_gust = Beta_bar_2_gust*(U(i)/b);
    Beta_3_gust = Beta_bar_3_gust*(U(i)/b);
    
    M_bar = M-qD*P2;
    C_bar = C-qD*P1;
    K_bar = K-qD*P0;
    
    s = [zeros(ndof), eye(ndof), zeros(ndof), zeros(ndof), zeros(ndof,1), zeros(ndof,1), zeros(ndof,1);
         -inv(M_bar)*K_bar, -inv(M_bar)*C_bar, qD*inv(M_bar)*P3, qD*inv(M_bar)*P4, qD*inv(M_bar)*P3_gust, qD*inv(M_bar)*P4_gust, qD*inv(M_bar)*P5_gust;
         zeros(ndof), eye(ndof), -Beta_1*eye(ndof), zeros(ndof), zeros(ndof,1), zeros(ndof,1), zeros(ndof,1);
         zeros(ndof), eye(ndof), zeros(ndof), -Beta_2*eye(ndof), zeros(ndof,1), zeros(ndof,1), zeros(ndof,1);
         zeros(1,ndof), zeros(1,ndof), zeros(1,ndof), zeros(1,ndof), -Beta_1_gust, 0, 0;
         zeros(1,ndof), zeros(1,ndof), zeros(1,ndof), zeros(1,ndof), 0, -Beta_2_gust, 0;
         zeros(1,ndof), zeros(1,ndof), zeros(1,ndof), zeros(1,ndof), 0, 0, -Beta_3_gust];
    Lambda(:,i) = reshape(eig(s),[],1);
end

% Root Locus Plot
figure() 
for i = 1:11
    plot(real(Lambda(i,:)),imag(Lambda(i,:)),'*'), hold on
end
xlabel('\sigma'),ylabel('j \omega'),title('Root Locus')
legend('location','NorthWest'),grid on

% % Flutter speed
% U_flutter_rl = zeros(11,1);
% for i = 1:11
%    U_flutter_rl(i,1) = interp1(real(Lambda(i,:)),U,0); 
% end
% min(U_flutter_rl)

%% Find the flutter speed and frequency of 2 DOF system using U-g method
k_inverse = linspace(0.01,7,50);
[g,omega,u] = ug_solve(a,b,c,ndof,M,K,G_str,rho,k_inverse);

g_h(:,1) = g(1,1,:);
g_alpha(:,1) = g(2,1,:);

omega_h(:,1) = (omega(1,1,:))./(2*pi);
omega_alpha(:,1) = (omega(2,1,:))./(2*pi);

u_h(:,1) = u(1,1,:);
u_alpha(:,1) = u(2,1,:);

figure(), hold all
plot(u_h,g_h), plot(u_alpha,g_alpha)
axis([0,35,-2,1]), xlabel("U (m/s)"), ylabel("g")
legend("h","\alpha",'Location','NorthEast'), title('U-g Method')

figure(), hold all
plot(u_h,omega_h), plot(u_alpha,omega_alpha)
axis([0,35,0,15]), xlabel("U (m/s)"), ylabel("\omega (Hz)"), 
legend("h","\alpha",'Location','NorthEast'), title('U-g Method')

% Find U when g=0 and omega at corresponding U
flutter_speeds = u_alpha(abs(g_alpha)<0.1);
flutter_speed = flutter_speeds(end-1);
flutter_frequency = omega_alpha(u_alpha==flutter_speed);

%% Frequency dependent response for gust velocity amplitude of 1 WG
U_flutter = flutter_speed;
U_flutter_fraction = 0.5; %[0.5,0.6,0.7,0.8,0.9,0.95]; %Fraction of U_flutter
natural_frequencies = zeros(2,2,length(U_flutter_fraction));
omega_sweep = 0:0.001:200;

for fraction = 1:length(U_flutter_fraction)
    qD_1 = 0.5*rho*(U_flutter_fraction(fraction)*U_flutter)^2;
    for i = 1:length(omega_sweep)
        reduced_frequency = (omega_sweep(i)*b)/(U_flutter_fraction(fraction)*U_flutter);
        [F_gust_1, F_gust_modified_1, H_1, Phi_1] = Gust_vector(reduced_frequency,a,b);
        [FF_1,GG_1] = NewTheodorsenFunction(reduced_frequency);
        [T_1] = NewTfunctions(a,c);
        [AA_1(:,:)] = NewAeroMatrix2DOF(reduced_frequency,a,c,b,T_1,FF_1,GG_1,ndof);
        
        A_1 = (-omega_sweep(i)^2)*M+1i*omega_sweep(i)*C+K-qD_1*AA_1; 
        B_1 = qD_1*F_gust_1*(1/(U_flutter_fraction(fraction)*U_flutter)); % Gust input
        y_1 = A_1\B_1;
        mag_1(:,i) = abs(y_1);
    end
    for i = 1:2
        pks_1 = findpeaks(mag_1(i,:));
        for j = 1:2
            % store 2 peak indices (h in first row and alpha in second)
            pks_1_matrix(i,j) = find(pks_1(j) == mag_1(i,:));
        end
    end

    string = ["h", "\alpha"];
    figure()
    for i = 1:2
        subplot(2,1,i)
        plot(omega_sweep, mag_1(i,:))
        xlabel('\omega (rad/s)')
        ylabel(strcat('Amplitude (', string(i),' / F)'))
        title(strcat(string(i)," ", 'response at', " ", num2str(U_flutter_fraction(fraction)), '*U_{flutter}'))
    end
    
    omega_1 = [omega_sweep(pks_1_matrix(1,1)), omega_sweep(pks_1_matrix(1,2));
               omega_sweep(pks_1_matrix(2,1)), omega_sweep(pks_1_matrix(2,2))];

    natural_frequencies(:,:,fraction) = omega_1;
end

% Plot of resonant frequencies as functions of flight speeds
figure()
plot([0.5,0.6,0.7,0.8,0.9,0.95],[29.105,29.884,30.966,32.530,34.969,36.661])
xlabel('U_{flutter}'),ylabel('Resonant Frequency (rad/s)')
title('Resonant Frequency (h) vs Flight Speed')

figure()
plot([0.5,0.6,0.7,0.8,0.9,0.95],[29.400,30.287,31.485,33.143,35.509,36.940])
xlabel('U_{flutter}'),ylabel('Resonant Frequency (rad/s)')
title('Resonant Frequency (\alpha) vs Flight Speed')

%% Plot response to a unit step vertical gust velocity
U_infinity = [0.7,0.8,0.9];
for i = 1:length(U_infinity)
    U_inf = U_infinity(i)*U_flutter;
    qD_inf = 0.5*rho*U_inf^2;
    
    P0 = P0_bar;
    P1 = P1_bar*(b/U_inf);
    P2 = P2_bar*(b/U_inf)^2;
    P3 = P3_bar;
    P4 = P4_bar;
    Beta_1 = Beta_bar_1*(U_inf/b);
    Beta_2 = Beta_bar_2*(U_inf/b);
    
    P0_gust = P0_bar_gust;
    P1_gust = P1_bar_gust*(b/U_inf);
    P3_gust = P3_bar_gust;
    P4_gust = P4_bar_gust;
    P5_gust = P5_bar_gust;
    Beta_1_gust = Beta_bar_1_gust*(U_inf/b);
    Beta_2_gust = Beta_bar_2_gust*(U_inf/b);
    Beta_3_gust = Beta_bar_3_gust*(U_inf/b);
    
    M_bar = M-(qD_inf)*P2;
    C_bar = C-(qD_inf)*P1;
    K_bar = K-(qD_inf)*P0;
    
    A_matrix = [zeros(ndof), eye(ndof), zeros(ndof), zeros(ndof), zeros(ndof,1), zeros(ndof,1), zeros(ndof,1);
                -inv(M_bar)*K_bar, -inv(M_bar)*C_bar, qD_inf*inv(M_bar)*P3, qD_inf*inv(M_bar)*P4, qD_inf*inv(M_bar)*P3_gust, qD_inf*inv(M_bar)*P4_gust, qD_inf*inv(M_bar)*P5_gust;
                zeros(ndof), eye(ndof), -Beta_1*eye(ndof), zeros(ndof), zeros(ndof,1), zeros(ndof,1), zeros(ndof,1);
                zeros(ndof), eye(ndof), zeros(ndof), -Beta_2*eye(ndof), zeros(ndof,1), zeros(ndof,1), zeros(ndof,1);
                zeros(1,ndof), zeros(1,ndof), zeros(1,ndof), zeros(1,ndof), -Beta_1_gust, 0, 0;
                zeros(1,ndof), zeros(1,ndof), zeros(1,ndof), zeros(1,ndof), 0, -Beta_2_gust, 0;
                zeros(1,ndof), zeros(1,ndof), zeros(1,ndof), zeros(1,ndof), 0, 0, -Beta_3_gust];
    B_matrix = [0,0;
                0,0;
                (qD_inf/U_inf)*inv(M_bar)*P0_gust,(qD_inf/U_inf)*inv(M_bar)*P1_gust;
                0,0;
                0,0;
                0,0;
                0,0;
                0,1/U_inf;
                0,1/U_inf;
                0,1/U_inf];
    C_matrix = [1,0,0,0,0,0,0,0,0,0,0;
                0,1,0,0,0,0,0,0,0,0,0];
    D_matrix = zeros(size(C_matrix,1),size(B_matrix,2));

    sys = ss(A_matrix,B_matrix,C_matrix,D_matrix);
    figure()
    step(sys,5)
    title(strcat('Response to Unit Step Gust at'," ",num2str(U_infinity(i)),'*U_{flutter}'))
end

%% U-g and F_gust Functions (other functions are in folder)
function [g,omega,u] = ug_solve(a,b,c,ndof,M,K,G_str,rho,k_inverse)
    T = NewTfunctions(a,c);
    K_bar = (eye(ndof) + 1i.*G_str)*K;
    B_bar = -K_bar;
    
    k = 1./k_inverse;
    for i = 1:length(k)
        [FF(i),GG(i)] = NewTheodorsenFunction(k(i));
        AA(:,:,i) = NewAeroMatrix2DOF(k(i),a,c,b,T,FF(i),GG(i),ndof);
        AA_bar(:,:,i) = -1*M-((0.5*rho*b^2)./(k(i).^2)).*AA(:,:,i);
        [V(:,:,i),D(:,:,i)] = eig(AA_bar(:,:,i)/(B_bar));
        Lambda(:,:,i) = D(:,:,i)*ones(ndof,1);
        g(:,:,i) = (imag(Lambda(:,:,i))./real(Lambda(:,:,i)));
        omega(:,:,i) = abs(sqrt(1./real(Lambda(:,:,i))));
        u(:,:,i) = (omega(:,:,i).*b)./k(i);
    end
end

function [F_gust, F_gust_modified, H, Phi] = Gust_vector(k,a,b)
    epsk=0.0000000001; % tolerance below which k is considered zero and C(k)=1
    for i = 1:length(k)
        % Bessel's Function 
        if k(i) <= epsk 
            J0=besselj(0,epsk);
            J1=besselj(1,epsk);
            Y0=bessely(0,epsk);
            Y1=bessely(1,epsk);  
        else
            J0=besselj(0,k(i));
            J1=besselj(1,k(i));
            Y0=bessely(0,k(i));
            Y1=bessely(1,k(i));
        end 
        % Theodorsen Function 
        [FF(i),GG(i)] = NewTheodorsenFunction(k(i));
        Ck(i) = FF(i)+1i*GG(i);
        % Sears Function 
        Phi(i) = (J0-1i*J1)*Ck(i)+1i*J1;
        % Modified Sears Function 
        H(i) = Phi(i)*exp(-1i*k(i));
    end
    % Gust Vector
    F_gust = 2*pi*diag([2*b (2*b)^2])*[-1; 1/2*(1/2+a)]*Phi;
    F_gust_modified = 2*pi*diag([2*b (2*b)^2])*[-1; 1/2*(1/2+a)]*H;
end