% 2D Wing+Control Surface
clear all; close all; clc

% Input Parameters
m_h        = 3.3930;
m_alpha    = 0.0135;
m_Beta     = 0.0003;
zeta_h     = 0.0113;
zeta_alpha = 0.01626;
zeta_Beta  = 0.0115;
K_h        = 2818.4;
K_alpha    = 37.34; 
K_Beta     = 3.895;
M          = [3.393, 0.0859, 0.0040;
             0.0859, 0.0135, 0.0008;
             0.0040, 0.0008, 0.0003];
K          = [2818.4, 0, 0;
             0, 37.34, 0;
             0, 0, 3.895];

% Natural Frequencies
omega_h             = sqrt(K_h/m_h)
omega_alpha         = sqrt(K_alpha/m_alpha)
omega_Beta          = sqrt(K_Beta/m_Beta)

% Damping Matrix
C_h     = 2*zeta_h*omega_h*m_h;
C_alpha = 2*zeta_alpha*omega_alpha*m_alpha;
C_Beta  = 2*zeta_Beta*omega_Beta*m_Beta;
C       = diag([C_h C_alpha C_Beta])

% Frequency Response Plots (Uncoupled)
omega = 0:0.001:200;

figure(1)
y1 = sqrt(((1./m_h).^2)./(((omega_h.^2-omega.^2).^2)+4.*(zeta_h.*omega_h.*omega).^2));
[max_y1, max_y1_index] = max(y1); % find max amplitude and natural frequency
omega_h_graph = omega(max_y1_index) % print natural frequency
idx1 = find(abs(y1-max(y1)/sqrt(2))<0.000005); % frequency at 1/sqrt(2) max amplitude
omega_h_delta = omega(idx1(2))-omega(idx1(1)); % frequency delta
zeta_h_graph = omega_h_delta/(2*omega_h_graph) % print damping ratio
plot(omega, y1) % plot frequency response
title('Frequency Response (h)'), xlabel('\omega (rad/s)'), ylabel('Amplitude')

figure(2)
y2 = sqrt(((1./m_alpha).^2)./(((omega_alpha.^2-omega.^2).^2)+4.*(zeta_alpha.*omega_alpha.*omega).^2));
[max_y2, max_y2_index] = max(y2); % find max amplitude and natural frequency
omega_alpha_graph = omega(max_y2_index) % print natural frequency
idx2 = find(abs(y2-max(y2)/sqrt(2))<0.00015); % frequency at 1/sqrt(2) max amplitude
omega_alpha_delta = omega(idx2(2))-omega(idx2(1)); % frequency delta
zeta_alpha_graph = omega_alpha_delta/(2*omega_alpha_graph) % print damping ratio
plot(omega, y2) % plot frequency response
title('Frequency Response (\alpha)'), xlabel('\omega (rad/s)'), ylabel('Amplitude')

figure(3)
y3 = sqrt(((1./m_Beta).^2)./(((omega_Beta.^2-omega.^2).^2)+4.*(zeta_Beta.*omega_Beta.*omega).^2));
[max_y3, max_y3_index] = max(y3); % find max amplitude and natural frequency
omega_Beta_graph = omega(max_y3_index) % print natural frequency
idx3 = find(abs(y3-max(y3)/sqrt(2))<0.0012); % frequency at 1/sqrt(2) max amplitude
omega_Beta_delta = omega(idx3(2))-omega(idx3(1)); % frequency delta
zeta_Beta_graph = omega_Beta_delta/(2*omega_Beta_graph) % print damping ratio
plot(omega, y3) % plot frequency response
title('Frequency Response (\beta)'), xlabel('\omega (rad/s)'), ylabel('Amplitude')

% Time Histories
t = [0,5]; % initial and end times

figure(4)
y1_0 = [1;0]; % initial perturbation
f1 = @(t,y)[y(2); (-K_h/m_h)*y(1)-(C_h/m_h)*y(2)];
[ts1,ys1] = ode45(f1,t,y1_0);
plot(ts1,ys1(:,1))
title('Time History (h)'), xlabel('Time (s)'), ylabel('Amplitude')
pks1 = findpeaks(ys1(:,1)); % find peaks
zeta_h_decrement = log(pks1(4)/pks1(5))/(2*pi) % used 4th and 5th peak

figure(5)
y2_0 = [1;0]; % initial perturbation
f2 = @(t,y)[y(2); (-K_alpha/m_alpha)*y(1)-(C_alpha/m_alpha)*y(2)];
[ts2,ys2] = ode45(f2,t,y2_0);
plot(ts2,ys2(:,1))
title('Time History (\alpha)'), xlabel('Time (s)'), ylabel('Amplitude')
pks2 = findpeaks(ys2(:,1)); % find peaks
zeta_alpha_decrement = log(pks2(4)/pks2(5))/(2*pi) % used 4th and 5th peak

figure(6)
y3_0 = [1;0]; % initial perturbation
f3 = @(t,y)[y(2); (-K_Beta/m_Beta)*y(1)-(C_Beta/m_Beta)*y(2)];
[ts3,ys3] = ode45(f3,t,y3_0);
plot(ts3,ys3(:,1))
title('Time History (\beta)'), xlabel('Time (s)'), ylabel('Amplitude')
pks3 = findpeaks(ys3(:,1)); % find peaks
zeta_Beta_decrement = log(pks3(4)/pks3(5))/(2*pi) % used 4th and 5th peak

% Coupled Response
% Input of [1; 0; 0]
string = ["h", "\alpha", "\beta"];
for i = 1:length(omega)
    A_1 = (-omega(i)^2)*M+1i*omega(i)*C+K; 
    b_1 = [1; 0; 0];
    y_1 = A_1\b_1;
    mag_1(:,i) = abs(y_1);
end
for i = 1:3
    pks_1 = findpeaks(mag_1(i,:));
    for k = 1:3
        % store 3 peak indices (h in first row, alpha in second, etc.)
        pks_1_matrix(i,k) = find(pks_1(k) == mag_1(i,:));
    end
end
figure(7)
for i = 1:3
    subplot(3,1,i)
    plot(omega, mag_1(i,:))
    xlabel('\omega (rad/s)')
    ylabel(strcat("Amplitude (", string(i)," / F)"))
end
omega_1 = [omega(pks_1_matrix(1,1)), omega(pks_1_matrix(1,2)), omega(pks_1_matrix(1,3));
           omega(pks_1_matrix(2,1)), omega(pks_1_matrix(2,2)), omega(pks_1_matrix(2,3));
           omega(pks_1_matrix(3,1)), omega(pks_1_matrix(3,2)), omega(pks_1_matrix(3,3))]

pks_1_1 = findpeaks(mag_1(1,:));
index_h_1_1 = find(abs(mag_1(1,1:40000)-pks_1_1(1,1)/sqrt(2))<0.00001); % frequency at 1/sqrt(2) max amplitude
omega_h_1_1_delta = omega(index_h_1_1(2))-omega(index_h_1_1(1)); % frequency delta
zeta_h_1_1_graph = omega_h_1_1_delta/(2*omega(pks_1_matrix(1,1))); % find damping ratio
index_h_2_1 = find(abs(mag_1(1,40000:90000)-pks_1_1(1,2)/sqrt(2))<0.0000001); % frequency at 1/sqrt(2) max amplitude
omega_h_2_1_delta = omega(index_h_2_1(2))-omega(index_h_2_1(1)); % frequency delta
zeta_h_2_1_graph = omega_h_2_1_delta/(2*omega(pks_1_matrix(1,2))); % find damping ratio
index_h_3_1 = find(abs(mag_1(1,120000:140000)-pks_1_1(1,3)/sqrt(2))<0.0000015625); % frequency at 1/sqrt(2) max amplitude
omega_h_3_1_delta = omega(index_h_3_1(2))-omega(index_h_3_1(1)); % frequency delta
zeta_h_3_1_graph = omega_h_3_1_delta/(2*omega(pks_1_matrix(1,3))); % find damping ratio

pks_2_1 = findpeaks(mag_1(2,:));
index_alpha_1_1 = find(abs(mag_1(2,1:40000)-pks_2_1(1,1)/sqrt(2))<0.000025); % frequency at 1/sqrt(2) max amplitude
omega_alpha_1_1_delta = omega(index_alpha_1_1(2))-omega(index_alpha_1_1(1)); % frequency delta
zeta_alpha_1_1_graph = omega_alpha_1_1_delta/(2*omega(pks_1_matrix(2,1))); % find damping ratio
index_alpha_2_1 = find(abs(mag_1(2,40000:90000)-pks_2_1(1,2)/sqrt(2))<0.0000045); % frequency at 1/sqrt(2) max amplitude
omega_alpha_2_1_delta = omega(index_alpha_2_1(2))-omega(index_alpha_2_1(1)); % frequency delta
zeta_alpha_2_1_graph = omega_alpha_2_1_delta/(2*omega(pks_1_matrix(2,2))); % find damping ratio
index_alpha_3_1 = find(abs(mag_1(2,120000:140000)-pks_2_1(1,3)/sqrt(2))<0.0000005); % frequency at 1/sqrt(2) max amplitude
omega_alpha_3_1_delta = omega(index_alpha_3_1(7))-omega(index_alpha_3_1(6)); % frequency delta
zeta_alpha_3_1_graph = omega_alpha_3_1_delta/(2*omega(pks_1_matrix(2,3))); % find damping ratio

pks_3_1 = findpeaks(mag_1(3,:));
index_Beta_1_1 = find(abs(mag_1(3,1:40000)-pks_3_1(1,1)/sqrt(2))<0.00001); % frequency at 1/sqrt(2) max amplitude
omega_Beta_1_1_delta = omega(index_Beta_1_1(2))-omega(index_Beta_1_1(1)); % frequency delta
zeta_Beta_1_1_graph = omega_Beta_1_1_delta/(2*omega(pks_1_matrix(3,1))); % find damping ratio
index_Beta_2_1 = find(abs(mag_1(3,40000:90000)-pks_3_1(1,2)/sqrt(2))<0.000003); % frequency at 1/sqrt(2) max amplitude
omega_Beta_2_1_delta = omega(index_Beta_2_1(2))-omega(index_Beta_2_1(1)); % frequency delta
zeta_Beta_2_1_graph = omega_Beta_2_1_delta/(2*omega(pks_1_matrix(3,2))); % find damping ratio
index_Beta_3_1 = find(abs(mag_1(3,120000:140000)-pks_3_1(1,3)/sqrt(2))<0.00000065); % frequency at 1/sqrt(2) max amplitude
omega_Beta_3_1_delta = omega(index_Beta_3_1(2))-omega(index_Beta_3_1(1)); % frequency delta
zeta_Beta_3_1_graph = omega_Beta_3_1_delta/(2*omega(pks_1_matrix(3,3))); % find damping ratio

zeta_1 = [zeta_h_1_1_graph, zeta_h_2_1_graph, zeta_h_3_1_graph;
          zeta_alpha_1_1_graph, zeta_alpha_2_1_graph, zeta_alpha_3_1_graph;
          zeta_Beta_1_1_graph, zeta_Beta_2_1_graph, zeta_Beta_3_1_graph]

% Input of [0; 1; 0]
for i = 1:length(omega)
    A_2 = (-omega(i)^2)*M+1i*omega(i)*C+K; 
    b_2 = [0; 1; 0];
    y_2 = A_2\b_2;
    mag_2(:,i) = abs(y_2);
end
for i = 1:3
    pks_2 = findpeaks(mag_2(i,:));
    for k = 1:3
        % store 3 peak indices (h in first row, alpha in second, etc.)
        pks_2_matrix(i,k) = find(pks_2(k) == mag_2(i,:));
    end
end
figure(8)
for i = 1:3
    subplot(3,1,i)
    plot(omega, mag_2(i,:))
    xlabel('\omega (rad/s)')
    ylabel(strcat("Amplitude (", string(i)," / F)"))
end
omega_2 = [omega(pks_2_matrix(1,1)), omega(pks_2_matrix(1,2)), omega(pks_2_matrix(1,3));
           omega(pks_2_matrix(2,1)), omega(pks_2_matrix(2,2)), omega(pks_2_matrix(2,3));
           omega(pks_2_matrix(3,1)), omega(pks_2_matrix(3,2)), omega(pks_2_matrix(3,3))]

pks_1_2 = findpeaks(mag_2(1,:));
index_h_1_2 = find(abs(mag_2(1,1:40000)-pks_1_2(1,1)/sqrt(2))<0.00002); % frequency at 1/sqrt(2) max amplitude
omega_h_1_2_delta = omega(index_h_1_2(2))-omega(index_h_1_2(1)); % frequency delta
zeta_h_1_2_graph = omega_h_1_2_delta/(2*omega(pks_2_matrix(1,1))); % find damping ratio
index_h_2_2 = find(abs(mag_2(1,40000:90000)-pks_1_2(1,2)/sqrt(2))<0.000004); % frequency at 1/sqrt(2) max amplitude
omega_h_2_2_delta = omega(index_h_2_2(2))-omega(index_h_2_2(1)); % frequency delta
zeta_h_2_2_graph = omega_h_2_2_delta/(2*omega(pks_2_matrix(1,2))); % find damping ratio
index_h_3_2 = find(abs(mag_2(1,120000:140000)-pks_1_2(1,3)/sqrt(2))<0.0000005); % frequency at 1/sqrt(2) max amplitude
omega_h_3_2_delta = omega(index_h_3_2(7))-omega(index_h_3_2(6)); % frequency delta
zeta_h_3_2_graph = omega_h_3_2_delta/(2*omega(pks_2_matrix(1,3))); % find damping ratio

pks_2_2 = findpeaks(mag_2(2,:));
index_alpha_1_2 = find(abs(mag_2(2,1:40000)-pks_2_2(1,1)/sqrt(2))<0.0001); % frequency at 1/sqrt(2) max amplitude
omega_alpha_1_2_delta = omega(index_alpha_1_2(4))-omega(index_alpha_1_2(3)); % frequency delta
zeta_alpha_1_2_graph = omega_alpha_1_2_delta/(2*omega(pks_2_matrix(2,1))); % find damping ratio
index_alpha_2_2 = find(abs(mag_2(2,40000:90000)-pks_2_2(1,2)/sqrt(2))<0.0001); % frequency at 1/sqrt(2) max amplitude
omega_alpha_2_2_delta = omega(index_alpha_2_2(2))-omega(index_alpha_2_2(1)); % frequency delta
zeta_alpha_2_2_graph = omega_alpha_2_2_delta/(2*omega(pks_2_matrix(2,2))); % find damping ratio
index_alpha_3_2 = find(abs(mag_2(2,120000:140000)-pks_2_2(1,3)/sqrt(2))<0.000005); % frequency at 1/sqrt(2) max amplitude
omega_alpha_3_2_delta = omega(index_alpha_3_2(2))-omega(index_alpha_3_2(1)); % frequency delta
zeta_alpha_3_2_graph = omega_alpha_3_2_delta/(2*omega(pks_2_matrix(2,3))); % find damping ratio

pks_3_2 = findpeaks(mag_2(3,:));
index_Beta_1_2 = find(abs(mag_2(3,1:40000)-pks_3_2(1,1)/sqrt(2))<0.00003); % frequency at 1/sqrt(2) max amplitude
omega_Beta_1_2_delta = omega(index_Beta_1_2(2))-omega(index_Beta_1_2(1)); % frequency delta
zeta_Beta_1_2_graph = omega_Beta_1_2_delta/(2*omega(pks_2_matrix(3,1))); % find damping ratio
index_Beta_2_2 = find(abs(mag_2(3,40000:90000)-pks_3_2(1,2)/sqrt(2))<0.0001); % frequency at 1/sqrt(2) max amplitude
omega_Beta_2_2_delta = omega(index_Beta_2_2(2))-omega(index_Beta_2_2(1)); % frequency delta
zeta_Beta_2_2_graph = omega_Beta_2_2_delta/(2*omega(pks_2_matrix(3,2))); % find damping ratio
index_Beta_3_2 = find(abs(mag_2(3,120000:140000)-pks_3_2(1,3)/sqrt(2))<0.00005); % frequency at 1/sqrt(2) max amplitude
omega_Beta_3_2_delta = omega(index_Beta_3_2(2))-omega(index_Beta_3_2(1)); % frequency delta
zeta_Beta_3_2_graph = omega_Beta_3_2_delta/(2*omega(pks_2_matrix(3,3))); % find damping ratio

zeta_2 = [zeta_h_1_2_graph, zeta_h_2_2_graph, zeta_h_3_2_graph;
          zeta_alpha_1_2_graph, zeta_alpha_2_2_graph, zeta_alpha_3_2_graph;
          zeta_Beta_1_2_graph, zeta_Beta_2_2_graph, zeta_Beta_3_2_graph]
      