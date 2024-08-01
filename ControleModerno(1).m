clear 
clc
close all

% Definir os parâmetros do sistema
k_p = 105e3;  % N/m
k_m = 30.25e3;  % N/m
k_cinto = 52.5e3;  % N/m
K = 21.678;  % Nm/rad
k_b = 40e3;  % N/m
b_b = 220; % Ns/m
b_p = 325;  % Ns/m
b_m = 2408.04;  % Ns/m
b_cinto = 250;  % Ns/m
B_rot = 2.4692;  % Nms/rad
l = 0.19;  % m
M_total = 70;  % kg
M_t = 42.609;  % kg
M_m = 21.238;  % kg
M_c = 6.153;  % kg
J_c = 0.02;  % kgm^2
M_b = 35;  % kg
g = 9.81; % m/s^2

% Definindo as matrizes
n1 = ((1 + M_c * l^2) / (J_c + M_c * l^2) / (M_c + M_m))^(-1);
n2 = (-M_c * l / (J_c + M_c * l^2));
n3 = (M_c + M_m + ((M_c^2) * (l^2)) / (J_c + M_c * l^2));
A = [
    0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 0, 1;
    -(k_b + k_cinto) / M_b, k_cinto / M_b, 0, 0, -(b_b + b_cinto) / M_b, b_cinto / M_b, 0, 0;
    k_cinto * n3, -(k_cinto + k_m) * n3, k_m * n3, M_c * l * (M_c * l * g - K) / (J_c + M_c * l^2) * n3, b_cinto / (M_t + M_c), -(b_cinto + b_m) / (M_t + M_c), b_m / (M_t + M_c), -B_rot * n3 / (J_c + M_c * l^2);
    0, k_m / M_m, -(k_m + k_p) / M_m, 0, 0, b_m / M_m, -(b_m + b_p) / M_m, 0;
    n1 * n2 * k_cinto, n1 * n2 * (-k_cinto - k_m), k_m * n1 * n2, (M_c * l * g - K) / (J_c + M_c * l^2), n1 * n2 * b_cinto, (-b_cinto - b_m) * n1 * n2, b_m * n1 * n2, -B_rot / (J_c + M_c * l^2)
];

B = [
    0;
    0;
    0;
    0;
    1 / M_b;
    0;
    0;
    0
];

E = [
    0, 0;
    0, 0;
    0, 0;
    0, 0;
    k_b / M_b, b_b / M_b;
    0, 0;
    k_p / M_m, b_p / M_m;
    0, 0
];

C =  [
    1, 0, 0, 0, 0, 0, 0, 0;
    -(k_b + k_cinto) / M_b, k_cinto / M_b, 0, 0, -(b_b + b_cinto) / M_b, b_cinto / M_b, 0, 0;
   ];

D =  [
    0;
    1 / M_b;
];

%Definindo o espaço de estados
sistema  = ss(A,B,C,D);


%Alocação de polos
P_alocacao = 1e3*[
  -0.0028 + 0.0000i;
  -0.0146 + 0.0453i;
  -0.0146 - 0.0453i;
  -0.0211 + 0.0818i;
  -0.0211 - 0.0818i;
  -0.0294 + 0.0000i;
  -0.0644 + 5.1173i;
  -0.0644 - 5.1173i       
];

K_alocacao = place(A, B , P_alocacao);

%LQR
Q = [
     1000,0,0,0,0,0,0,0;
     0,1500,0,0,0,0,0,0;
     0,0,2000,0,0,0,0,0;
     0,0,0,5000,0,0,0,0;
     0,0,0,0,7000,0,0,0;
     0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,500,0;
     0,0,0,0,0,0,0,0
]; 

R = 0.01; 

[K_lqr,S_lqr,P_lqr] = lqr(A, B, Q, R);

%Simulação
tempoDeSimulacao = 2; %s
amostra = 0.001;
tempo = 0:amostra:tempoDeSimulacao;
x0 = [0; 0; 0; 0; 0; 0; 0; 0];

[t,x] = ode45(@(t,x) SemControle(t,x,A,B,E), tempo, x0);
[t_alocacao,x_alocacao] = ode45(@(t,x) AlocacaoDePolos(t,x,A,B,K_alocacao,E), tempo, x0);
[t_lqr,x_lqr] = ode45(@(t,x) LQR(t,x,A,B,K_lqr,E), tempo, x0);
EntradaDeControle_alocacao = K_alocacao * x';
EntradaDeControle_lqr = K_lqr * x';

% Extraindo as velocidades
%Sem controle
xb_dot = x(:, 5); 
xt_dot = x(:, 6); 
xm_dot = x(:, 7); 
theta_dot = x(:, 8);
%Alocação de polos
xb_dot_alocacao = x_alocacao(:, 5); 
xt_dot_alocacao = x_alocacao(:, 6); 
xm_dot_alocacao = x_alocacao(:, 7); 
theta_dot_alocacao = x_alocacao(:, 8);
%LQR
xb_dot_lqr = x_lqr(:, 5); 
xt_dot_lqr = x_lqr(:, 6); 
xm_dot_lqr = x_lqr(:, 7); 
theta_dot_lqr = x_lqr(:, 8);

% Calculando as acelerações através das derivadas das velocidades
%Sem controle
xb_ddot = gradient(xb_dot, t);
xt_ddot = gradient(xt_dot, t);
xm_ddot = gradient(xm_dot, t);
theta_ddot = gradient(theta_dot, t);
%Alocação de polos
xb_ddot_alocacao = gradient(xb_dot_alocacao, t_alocacao);
xt_ddot_alocacao = gradient(xt_dot_alocacao, t_alocacao);
xm_ddot_alocacao = gradient(xm_dot_alocacao, t_alocacao);
theta_ddot_alocacao = gradient(theta_dot_alocacao, t_alocacao);
%LQR
xb_ddot_lqr = gradient(xb_dot_lqr, t_lqr);
xt_ddot_lqr = gradient(xt_dot_lqr, t_lqr);
xm_ddot_lqr = gradient(xm_dot_lqr, t_lqr);
theta_ddot_lqr = gradient(theta_dot_lqr, t_lqr);

% Acelerações calculadas
%Sem controle
aceleracoes = [xb_ddot, xt_ddot, xm_ddot, theta_ddot]/9.81;
%Alocação de polos
aceleracoes_alocacao = [xb_ddot_alocacao, xt_ddot_alocacao, xm_ddot_alocacao, theta_ddot_alocacao]/9.81;
%LQR
aceleracoes_lqr = [xb_ddot_lqr, xt_ddot_lqr, xm_ddot_lqr, theta_ddot_lqr]/9.81;

%Simulando
figure(1)
plot(t, x(:, 1), t_alocacao, x_alocacao(:,1), t_lqr, x_lqr(:,1));
xlabel('Tempo (s)');
ylabel('Deslocamento (m)');
legend('Sem controle','Alocação de Polos','LQR')
title('Deslocamento do Banco');


figure(2)
plot(t, aceleracoes(:, 1), t, aceleracoes_alocacao(:, 1), t, aceleracoes_lqr(:, 1));
xlabel('Tempo (s)');
ylabel('Aceleração (Gs)');
title('Aceleraçao do Banco');

line([0.04, 0.1], [-37, -37], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
xline(0.04, 'LineWidth', 2);
xline(0.1, 'LineWidth', 2);
ylim([-44 -20])
xlim([0 0.15])

legend('Sem controle','Alocação de Polos','LQR','Location','northeast','Orientation','vertical')

figure(3)
plot(t,EntradaDeControle_alocacao, t, EntradaDeControle_lqr)
legend('Alocação de Polos', 'LQR')
xlabel('Tempo(s)')
ylabel('u(t)')
title('Entrada de controle')

figure(4)
plot(t, aceleracoes(:, 1));
xlabel('Tempo (s)');
ylabel('Aceleração (Gs)');
title('Aceleraçao do Banco');

figure(5)
plot(t_alocacao, x_alocacao(:,1));
xlabel('Tempo (s)');
ylabel('Deslocamento (m)');
title('Deslocamento do Banco por alocação de polos');

figure(6)
plot(t_alocacao, aceleracoes_alocacao(:, 1));
xlabel('Tempo (s)');
ylabel('Aceleração (Gs)');
title('Aceleraçao do Banco por alocação de polos');

figure(7)
plot(t_lqr, x_lqr(:,1));
xlabel('Tempo (s)');
ylabel('Deslocamento (m)');
title('Deslocamento do Banco por LQR');

figure(8)
plot(t_alocacao, aceleracoes_lqr(:, 1));
xlabel('Tempo (s)');
ylabel('Aceleração (Gs)');
title('Aceleraçao do Banco por LQR');

figure(9)
plot(t, aceleracoes(:, 1), t_alocacao, aceleracoes_alocacao(:, 1), t_lqr, aceleracoes_lqr(:, 1));
xlabel('Tempo (s)');
ylabel('Aceleração (Gs)');
title('Aceleraçao do Banco');
legend('Sem controle', 'Alocação de Polos', 'LQR')

sys = zpk([], P_alocacao, 1);
figure(10);
pzmap(sys);
title('Diagrama de Polos do Sistema em Malha Fechada');
grid on;

MC = ctrb(A,B);
n = rank(MC, 0.00001);

if n == length(A)
    disp('Controlável')
else
    disp('Não controlável')
end

MO = obsv(A,C);
n = rank(MO, 0.00001);

if n == length(A)
    disp('Observável')
else
    disp('Não observável')
end

function f_t = velocidadeNaColisao(t)
    v = 50 / 3.6; % m/s
    p2 = v^2;
    p1 = (-p2) / 0.3;

    if isscalar(t)
        if t <= 0.3
            quadratic_part = p1 * t + p2;
            f_t = sqrt(max(quadratic_part, 0));
        else
            f_t = 0;
        end
    else
        f_t = zeros(size(t));
        quadratic_part = p1 * t + p2;
        f_t(t <= 0.3) = sqrt(max(quadratic_part(t <= 0.3), 0));
        f_t(t > 0.3) = 0;
    end
end
function dxdt = SemControle(t,x,A,B,E)
    r2 = velocidadeNaColisao(t); %m/s

    L = 0.38;
    k = 20;
    f = 0.01;
    r1 = @(t) L ./ (1 + exp(-k * (t - f))); %m

    c = @(t) 3 + sin(2*pi*t); %Ns/m
    
    u = c(t) .* r2; 
   
    w = [r1(t);
     r2
    ];
    
    dxdt = A * x + B* u + E * w;
end
function dxdt = AlocacaoDePolos(t,x,A,B,K_alocacao,E)
    ref = [
           0;
           0;
           0;
           0;
           0;
           0;
           0;
           0   
    ];
    r2 = velocidadeNaColisao(t); %m/s

    L = 0.38;
    k = 20;
    f = 0.01;
    r1 = @(t) L ./ (1 + exp(-k * (t - f))); %m

    w = [r1(t);
     r2
    ];

    dxdt = A * x + B * (-K_alocacao * (x - ref)) + E * w;
end
function dxdt = LQR(t,x,A,B,K_lqr,E)
    ref = [
           0;
           0;
           0;
           0;
           0;
           0;
           0;
           0   
    ];
    r2 = velocidadeNaColisao(t); %m/s

    L = 0.38;
    k = 20;
    f = 0.01;
    r1 = @(t) L ./ (1 + exp(-k * (t - f))); %m

    w = [r1(t);
     r2
    ];

    dxdt = A * x + B * (-K_lqr * (x - ref)) + E * w;
end