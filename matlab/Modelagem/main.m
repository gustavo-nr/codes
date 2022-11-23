clear all
clc
close all

%% Análise
% Para simular uma função, atribua um valor para a variável "code" conforme
% a enumeração listada abaixo

% 1. Matriz de funções de transferência
% 2. Diagramas de Bode
% 3. Simulações
% 4. Matriz de transição e termos forçantes - simulação e comparação com 
%    modelo linear (extra)

code = 4;

%% Dados do sistema
% Matrizes do espaço de estados
A = readmatrix("A_l.txt");
B = readmatrix("B_l.txt");
C = readmatrix("C_l.txt");
D = zeros(size(C,1),size(B,2));

%% Vetor de tempo
t_i = 0; % Tempo inicial
step = 0.001; % Passo de tempo
t_f = 2.5; % Tempo final
t = t_i:step:t_f;

%% Condições iniciais ( ex.: phi = phi_eq + delta_phi)
delta_phi = 0;
delta_psi = 0;
delta_theta = 0;
delta_phip = 0;
delta_psip = 0;
delta_thetap = 0;
delta_hx = 0;
delta_hy = 0;
delta_hz = 0;

%% Vetor de entradas
% - Para definir um degrau, basta inserir uma valor não nulo para uma das
% entradas e um intervalo de tempo significativo.

% - Para definir um impulso, basta inserir uma valor não nulo para uma
% das entradas e um intervalo de tempo infinitesimal.

mag= [0.1 0.1 0.1]; %Magnitude do sinal de entrada
ti = [0 0 0]; %Tempo incial 
tf = [0.5 0.5 0.5]; %Tempo final
u = u_selector(mag,ti,tf,t);

%% Ponto de operação
phi_eq = 0.01*pi/6;
psi_eq = 1*pi/6;
theta_eq = 1*pi/6;
phip_eq = 0;
psip_eq = 0;
thetap_eq = 0;
hx_eq = 0;
hy_eq = 0;
hz_eq = 0;

x_eq = [phi_eq psi_eq theta_eq phip_eq psip_eq thetap_eq hx_eq hy_eq hz_eq];
x0_lin = [delta_phi delta_psi delta_theta delta_phip delta_psip delta_thetap delta_hx delta_hy delta_hz];
x0_nlin = [phi_eq+delta_phi psi_eq+delta_psi theta_eq+delta_theta phip_eq+delta_phip psip_eq+delta_psip...
           thetap_eq+delta_thetap hx_eq+delta_hx hy_eq+delta_hy hz_eq+delta_hz];

%% 1. Matriz de funções de transferência
if code == 1
    % Para obter a matriz de função de transferência, basta analisar a
    % a variável simbólica TFM.
    
    n_round = 10;
    TFM = tf_matrix(A,B,C,D,n_round);
    TFM1_latex = latex(TFM(:,1));
    TFM2_latex = latex(TFM(:,2));
    TFM3_latex = latex(TFM(:,3));
end

%% 2. Diagramas de Bode
if code == 2
    % Para plotar os diagramas de Bode do sistema, basta definir 'code = 3'.
    
    bode_plot(A,B,C,D)
end

%% 3. Simulações
if code == 3
    % Os gráficos plotados são definidos na função sys_sim. Para qualquer
    % variação, sugere-se utilizar as saídas [y_lin,y_nlin,t] da função 
    % @sys_sim(A,B,C,D,x0_lin,x0_nlin,x_eq,u,t,mag,ti,tf).
    
    [y_lin,y_nlin,t] = sys_sim(A,B,C,D,x0_lin,x0_nlin,x_eq,u,t,mag,ti,tf);
end

%% 4. Matriz de transição e termos forçantes
if code == 4
    % Para obter as matrizes de transição (Phi) e de termos forçantes
    % (Gamma), basta analisar as variáveis [Phi_dt,Gamma_dt].
    
    [Phi_dt,Gamma_dt,x] = transition_matrix(A,B,t,x0_nlin',u');
    [y_lin,x,t] = sys_transition_sim(A,B,C,D,x0_lin,x0_nlin,x_eq,u,t);
end
