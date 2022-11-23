function [Phi_dt,Gamma_dt,x] = transition_matrix(A,B,t,x0,u)
%% Dados de entrada
    % Matriz de estados A
    dim_A = length(A);
    
    %% Vetor de tempo
    ti = t(1); % instante inicial
    tf = t(end); % instante final
    n = length(t);
    dt = (tf-ti)/(n-1); % passo
    
    %% Matriz de Transição
    Phi_dt = eye(dim_A) + A * dt + (A^2 * dt^2)/(2) + (A^3 * dt^3)/(6);
    
    %% Termo forçante
    Gamma_dt = dt * (eye(dim_A) + (A*dt)/2 + (A^2 * dt^2)/6);
    M_forc = Gamma_dt*B;
    
    %% Solução do sistema de equações
    x = zeros(dim_A,n);
    x(:,1) = x0;
    for i=1:n-1
        x(:,i+1) = Phi_dt * x(:,i) + M_forc*u(:,i);
    end
end