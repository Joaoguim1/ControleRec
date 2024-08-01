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
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    1 / M_b, k_b / M_b, b_b / M_b;
    0, 0, 0;
    0, k_p / M_m, b_p / M_m;
    0, 0, 0
];

C =  [
    1, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0;
    -(k_b + k_cinto) / M_b, k_cinto / M_b, 0, 0, -(b_b + b_cinto) / M_b, b_cinto / M_b, 0, 0;
    k_cinto * n3, -(k_cinto + k_m) * n3, k_m * n3, M_c * l * (M_c * l * g - K) / (J_c + M_c * l^2) * n3, b_cinto / (M_t + M_c), -(b_cinto + b_m) / (M_t + M_c), b_m / (M_t + M_c), -B_rot * n3 / (J_c + M_c * l^2);
    0, k_m / M_m, -(k_m + k_p) / M_m, 0, 0, b_m / M_m, -(b_m + b_p) / M_m, 0;
    n1 * n2 * k_cinto, n1 * n2 * (-k_cinto - k_m), k_m * n1 * n2, (M_c * l * g - K) / (J_c + M_c * l^2), n1 * n2 * b_cinto, (-b_cinto - b_m) * n1 * n2, b_m * n1 * n2, -B_rot / (J_c + M_c * l^2)
];

D =  [
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    1 / M_b, k_b / M_b, b_b / M_b;
    0, 0, 0;
    0, k_p / M_m, b_p / M_m;
    0, 0, 0
];

% Definindo o espaço de estados
EspacoDeEstados = ss(A, B, C, D);
FT = tf(EspacoDeEstados(5, 1));
G = FT;

% Traçar o lugar das raízes variando Kp
figure(1)
rlocus(G);  % Corrigido para traçar o lugar das raízes da planta G
grid on
title('Lugar das Raízes para Variação de Kp');

% Ganho proporcional ótimo
Kp = 0.55;

% Função de transferência em malha fechada com Kp
FT_MF = feedback(Kp * G, 1);

% Função de transferência integrativa
s = tf('s');
FT_I = FT_MF / s;

figure(2)
rlocus(FT_I);  % Traçar lugar das raízes do sistema com ação integrativa
grid on 
title('Lugar das Raízes para Variação de Ki');

% Ganho integral ótimo
Ki = 20;

% Função de transferência em malha fechada com Kp e Ki
FT_PI = feedback(Kp * G + Ki / s * G, 1);  % Ajustado para FT_PI

% Função de transferência derivativa
FT_d = s * FT_PI;

figure(3)
rlocus(FT_d);  % Traçar lugar das raízes do sistema com ação derivativa
grid on 
title('Lugar das Raízes para Variação de Kd');

% Ganho derivativo
Kd = 0;

% Controlador PID
PID = pid(Kp, Ki, 0);

% Função de transferência em malha fechada com o controlador PID
T = feedback(PID * G, 1);

figure(4)
margin(T);

tempoDeSimulacao = 2; %s
amostra = 0.0001;
tempo = 0:amostra:tempoDeSimulacao;
x0 = [0; 0; 0; 0; 0; 0; 0; 0];

u = 22432.86*ones(size(tempo));


[y, t] = lsim(G, u, tempo, x0);
[y_pi, t_pi] = lsim(T, u, tempo, x0);
y = y/9.81;
y_pi = y_pi/9.81;

figure(5)
plot(t, y, t_pi, y_pi)
xlabel('Tempo (s)');
ylabel('Aceleração (Gs)');
title('Aceleraçao do Banco');
line([0.04, 0.1], [-37, -37], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
xline(0.04, 'LineWidth', 2);
xline(0.1, 'LineWidth', 2);
ylim([-44 -20])
xlim([0 0.15])
legend('Malha aberta', 'PI')

[sysb,g] = balreal(FT);
disp(g);
sys_red = modred(sysb,[6 7 8], 'del');
tf_red = tf(sys_red);
[num_red, den_red] = tfdata(tf_red,'v');
disp(tf_red);

% Plotando o diagrama de Bode para as duas funções de transferência
figure(6);
margin(FT); % Diagrama de Bode da função de transferência original com margens
hold on;
margin(tf_red); % Diagrama de Bode da função de transferência reduzida com margens

title('Diagrama de Bode - Comparação entre Funções de Transferência');
grid on;

% Traçar o lugar das raízes variando Kp ordem reduzida
figure(7)
rlocus(tf_red);  % Corrigido para traçar o lugar das raízes da planta G
grid on
title('Lugar das Raízes para Variação de Kp red');

% Ganho proporcional ótimo reduzido
Kp_red = 0.47;

% Função de transferência em malha fechada com Kp reduzido
FT_MF_red = feedback(Kp_red * tf_red, 1);

% Função de transferência integrativa reduzida
s = tf('s');
FT_I_red = FT_MF_red / s;

figure(8)
rlocus(FT_I_red);  % Traçar lugar das raízes do sistema com ação integrativa reduzida
grid on 
title('Lugar das Raízes para Variação de Ki red');

% Ganho integral ótimo reduzido
Ki_red = 15.5;

% Função de transferência em malha fechada com Kp e Ki reduzidos
FT_PI_red = feedback(Kp_red * tf_red + Ki_red / s * tf_red, 1);  % Ajustado para FT_PI red

% Função de transferência derivativa reduzida
FT_d_red = s * FT_PI_red;

figure(9)
rlocus(FT_d_red);  % Traçar lugar das raízes do sistema com ação derivativa red
grid on 
title('Lugar das Raízes para Variação de Kd red');

% Ganho derivativo
Kd_red = 0;

% Controlador PID
PID_red = pid(Kp_red, Ki_red, Kd_red);

% Função de transferência em malha fechada com o controlador PID
T_red = feedback(PID_red * tf_red, 1);

figure(10)
margin(T_red);

[y_pi_red, t_pi_red] = lsim(T_red, u, tempo, x0);
y_pi_red = y_pi_red/9.81;

figure(11)
plot(t_pi_red, y_pi_red, t_pi, y_pi)
xlabel('Tempo (s)');
ylabel('Aceleração (Gs)');
title('Aceleraçao do Banco');
%line([0.04, 0.1], [-37, -37], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
%xline(0.04, 'LineWidth', 2);
%xline(0.1, 'LineWidth', 2);
%ylim([-44 -20])
%xlim([0 0.15])
legend('PI Ordem 5', 'PI Ordem 8')

% polos = pole(FT);
% zeros = zero(FT);
% format long
% 
% disp('Polos da função de transferência:');
% disp(polos);
% disp('Zeros da função de transferência:');
% disp(zeros);