%% Seguidor de referência constante

clc
clear
close

%% Importa matrizes

addpath('Matrizes\')
addpath("Imagens\Controle Moderno\")

A = importdata('matrix_A1lin.txt');
B = importdata('matrix_B1lin.txt');
C = importdata('matrix_C1.txt');
D = importdata('matrix_D1.txt');

%% Define a referência constante
x_ref = [pi/4 pi/6 pi/3 0 0 0];
seq = 'ZXZ';
q_ref = angle2quat(x_ref(1),x_ref(2),x_ref(3),seq);
xq_ref = [q_ref(2:4),x_ref(4:6)];

%% Regulador por alocação de polos
p1=-0.3+0.3i;
p2=-0.4+0.3i;
p = [p1 conj(p1) p2 conj(p2) -0.5 -0.4];
K = place(A,B,p);
F = A-B*K;

%% Cálculo da entrada em regime permanente
Gamma = [A,B;C,D];
B_ls = [zeros(6,6);eye(3,6);zeros(3,6)];
N = linsolve(Gamma,B_ls);
Nx = N(1:6,:);
Nu = N(7:9,:);
u_rp = (Nu+K*Nx)*xq_ref';
B2 = B;
B2(4:6,:) = B2(4:6,:).*u_rp;

%% Simulação do sistema
sys_aloc = ss(F,B2,C,D);
[y,t,x]=step(sys_aloc,30);

q0 = q_ref(1);
q1 = y(:,1,1);
q2 = y(:,2,2);
q3 = y(:,3,3);
w1 = y(:,4,1);
w2 = y(:,5,2);
w3 = y(:,6,3);

[psi,theta,phi] = quat2eulang(q0,q1,q2,q3,seq);

figure(1)
plot(t,w1,LineWidth=1.20)
hold on
plot(t,w2,LineWidth=1.20)
plot(t,w3,LineWidth=1.20)
grid on
title("Velocidade angular ao longo do tempo")
xlabel("Tempo [s]")
ylabel("Velocidade angular [rad/s]")
legend('w1','w2','w3')
hold off

baseFileName = sprintf('Image_%s.png', "seg_cte_w_t");
fullFileName = fullfile("Imagens\Controle Moderno\", baseFileName);
saveas(1, fullFileName);

figure(2)
plot(t,q1,LineWidth=1.20)
hold on
plot(t,q2,LineWidth=1.20)
plot(t,q3,LineWidth=1.20)
grid on
title("Posição angular em quatérnios ao longo do tempo")
xlabel("Tempo [s]")
ylabel("Posição angular [rad]")
legend('q1','q2','q3','Location','northwest')
hold off

baseFileName = sprintf('Image_%s.png', "seg_cte_quat_t");
fullFileName = fullfile("Imagens\Controle Moderno\", baseFileName);
saveas(2, fullFileName);

figure(3)
plot(t,psi,LineWidth=1.20)
hold on
plot(t,theta,LineWidth=1.20)
plot(t,phi,LineWidth=1.20)
grid on
title("Posição angular em ângulos de Euler ao longo do tempoo")
xlabel("Tempo [s]")
ylabel("Posição angular [rad]")
legend('psi','theta','phi','Location','southeast')
hold off

baseFileName = sprintf('Image_%s.png', "seg_cte_euler_t");
fullFileName = fullfile("Imagens\Controle Moderno", baseFileName);
saveas(3, fullFileName);



