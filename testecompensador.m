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

% Ganho proporcional e integral ótimos
Kp = 0.55;
Ki = 20;

% Função de transferência em malha fechada com Kp e Ki
s = tf('s');
FT_PI = feedback(Kp * G + Ki / s * G, 1);

% Coeficientes do numerador e denominador do sistema
Y = [0.02857, 5.526, 7.483e5, 4.914e7, 5.94e9, 1.522e11, 6.092e11, -0.09089, -2.809e-15];
X = [1, 206.8, 2.621e7, 2.044e9, 2.872e11, 1.052e13, 4.817e14, 9.119e15, 3.425e16];

% Definindo a função de transferência
H = tf(Y, X);

% Escolher um zero e um polo para o compensador
z_c = -5;  % Posição do zero do compensador
p_c = -9;  % Posição do polo do compensador

% Definir a função de transferência do compensador
Gc = tf([1, -z_c], [1, -p_c]);

% Sistema compensado
H_compensado = series(Gc, H);

% Calcular o ganho necessário para atingir os polos desejados
dominant_pole = -5.73 + 42.85i;  % exemplo de polo dominante do sistema original
eval_H_compensado = evalfr(H_compensado, dominant_pole);
K = 1 / abs(eval_H_compensado);

% Aplicar o ganho ao sistema compensado
H_compensado = K * H_compensado;

% Função de transferência em malha fechada do sistema compensado
T2 = feedback(H_compensado, 1);

% Simulação das respostas ao degrau dos dois sistemas
tempoDeSimulacao = 2; %s
amostra = 0.0001;
tempo = 0:amostra:tempoDeSimulacao;
x0 = [0; 0; 0; 0; 0; 0; 0; 0];
u = 22432.86*ones(size(tempo));

[y_pi, t_pi] = lsim(FT_PI, u, tempo, x0);
y_pi = y_pi/9.81;

[y_pi2, t_pi2] = lsim(T2, u, tempo, x0);
y_pi2 = y_pi2/9.81;

% Plot comparison of time responses
figure
plot(t_pi, y_pi, t_pi2, y_pi2)
ylim([-50 70])
xlim([0 0.8])
title('Comparação sistema em malha fechada com e sem compensador');
legend('PI - Lugar das raízes', 'PI com compensador')
ylabel('Aceleração (Gs)')
xlabel('Tempo(s)')