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
    0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0;
    -(k_b + k_cinto) / M_b, k_cinto / M_b, 0, 0, -(b_b + b_cinto) / M_b, b_cinto / M_b, 0, 0;
    k_cinto * n3, -(k_cinto + k_m) * n3, k_m * n3, M_c * l * (M_c * l * g - K) / (J_c + M_c * l^2) * n3, b_cinto / (M_t + M_c), -(b_cinto + b_m) / (M_t + M_c), b_m / (M_t + M_c), -B_rot * n3 / (J_c + M_c * l^2);
    0, k_m / M_m, -(k_m + k_p) / M_m, 0, 0, b_m / M_m, -(b_m + b_p) / M_m, 0;
    n1 * n2 * k_cinto, n1 * n2 * (-k_cinto - k_m), k_m * n1 * n2, (M_c * l * g - K) / (J_c + M_c * l^2), n1 * n2 * b_cinto, (-b_cinto - b_m) * n1 * n2, b_m * n1 * n2, -B_rot / (J_c + M_c * l^2)
];

D =  [
    0;
    0;
    0;
    0;
    1 / M_b;
    0;
    0;
    0
];

C1=[
    0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0;
    k_cinto * n3, -(k_cinto + k_m) * n3, k_m * n3, M_c * l * (M_c * l * g - K) / (J_c + M_c * l^2) * n3, b_cinto / (M_t + M_c), -(b_cinto + b_m) / (M_t + M_c), b_m / (M_t + M_c), -B_rot * n3 / (J_c + M_c * l^2);
    0, k_m / M_m, -(k_m + k_p) / M_m, 0, 0, b_m / M_m, -(b_m + b_p) / M_m, 0;
    n1 * n2 * k_cinto, n1 * n2 * (-k_cinto - k_m), k_m * n1 * n2, (M_c * l * g - K) / (J_c + M_c * l^2), n1 * n2 * b_cinto, (-b_cinto - b_m) * n1 * n2, b_m * n1 * n2, -B_rot / (J_c + M_c * l^2)

];

D1 = [
    0;
    0;
    0;
    0;
    0;
    0];
    

C2 = [
    1, 0, 0, 0, 0, 0, 0, 0;
    -(k_b + k_cinto) / M_b, k_cinto / M_b, 0, 0, -(b_b + b_cinto) / M_b, b_cinto / M_b, 0, 0;
];
D2 = [0;
     1 / M_b
];

C3 = [
    1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0;
    -(k_b + k_cinto) / M_b, k_cinto / M_b, 0, 0, -(b_b + b_cinto) / M_b, b_cinto / M_b, 0, 0;
    0, k_m / M_m, -(k_m + k_p) / M_m, 0, 0, b_m / M_m, -(b_m + b_p) / M_m, 0;
];
D3 = [0;
     0;
     1 / M_b;
     1 / M_m
];


C4=[
   1, 0, 0, 0, 0, 0, 0, 0;
   ];
D4 = [
    0;
    ]

%Definindo o espaço de estados
sistema  = ss(A,B,C,D);
FT1 = tf(sistema(5,1));
%% Realização Balanceada
[sys_bal, g] = balreal(sistema);

% Exibindo a realização balanceada e os gramianos
%disp('Realização Balanceada:');
%disp(sys_bal);
%disp('Gramianos:');
%disp(g);

%% Ordem Reduzida
% Definir um limiar para escolher os estados a serem mantidos
limiar = 0.1;
elim = find(g < limiar); % Índices dos estados a serem eliminados

% Reduzir o sistema
sys_red = modred(sys_bal, elim);

% Exibir o sistema reduzido
%disp('Sistema Reduzido:');
%disp(sys_red);


% Obtendo a função de transferência reduzida (FT1_red) a partir do sistema reduzido
FT1_red = tf(sys_red(5, 1));

% Exibir a função de transferência original e reduzida
%disp('Função de Transferência Original:');
%disp(FT1);

%disp('Função de Transferência Reduzida:');
%disp(FT1_red);

% Plotando o diagrama de Bode para as duas funções de transferência
%figure(9);
%bode(FT1, {1e-2, 1e3}); % Diagrama de Bode da função de transferência original
%hold on;
%bode(FT1_red, {1e-2, 1e3}); % Diagrama de Bode da função de transferência reduzida
%legend('Original', 'Reduzida');
%title('Diagrama de Bode - Comparação entre Funções de Transferência');
%grid on;

% Resposta ao impulso do sistema original
t = 0:0.001:2;  % Vetor de tempo de 0 a 2 segundos, com passo de 0.01 segundos

[y_original, t_original] = impulse(sistema, t);

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

%plotarpoloslqr
sistlqr = ss(A-B*K_lqr,B,C2,D2);
%figure(14);
%pzmap(sistlqr);
%title('Diagrama de Polos do Sistema em Malha Fechada por LQR');
%grid on;



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

% % Conversão de acelerações para G/
% aceleracoes_g = aceleracoes / 9.81;
% aceleracoes_alocacao_g = aceleracoes_alocacao / 9.81;
% aceleracoes_lqr_g = aceleracoes_lqr / 9.81;
% 
% % Identificação dos picos de aceleração e suas durações
% [picos, locs] = findpeaks(aceleracoes_g(:, 1), 'MinPeakProminence', 1);
% picos(end) = [];
% duracao_picos = diff(locs) * amostra; % duração em segundos
% 
% [picos_alocacao, locs_alocacao] = findpeaks(aceleracoes_alocacao_g(:, 1), 'MinPeakProminence', 1);
% picos_alocacao(end) = [];
% duracao_picos_alocacao = diff(locs_alocacao) * amostra;
% 
% [picos_lqr, locs_lqr] = findpeaks(aceleracoes_lqr_g(:, 1), 'MinPeakProminence', 1);
% picos_lqr(end) = [];
% duracao_picos_lqr = diff(locs_lqr) * amostra;
% 
% % Comparação com o gráfico de referência (simplificado)
% figure;
% hold on;
% plot(picos, duracao_picos, 'ro', 'DisplayName', 'Picos de Aceleração (Sem Controle)');
% plot(picos_alocacao, duracao_picos_alocacao, 'bo', 'DisplayName', 'Picos de Aceleração (Alocação de Polos)');
% plot(picos_lqr, duracao_picos_lqr, 'go', 'DisplayName', 'Picos de Aceleração (LQR)');
% xlabel('Aceleração de Pico (G)');
% ylabel('Duração do Pulso (s)');
% legend('show');
% title('Comparação de Danos ao Ser Humano');
% hold off;

%Simulando
%figure(1)
%plot(t, x(:, 1), t_alocacao, x_alocacao(:,1), t_lqr, x_lqr(:,1));
%xlim([0 .8])
%ylim([0.5 2])
%xlabel('Tempo (s)');
%ylabel('Deslocamento (m)');
%legend('Sem controle','Alocação de Polos','LQR')
%title('Deslocamento do Banco');


%figure(2)
%plot(t, aceleracoes(:, 1), t, aceleracoes_alocacao(:, 1), t, aceleracoes_lqr(:, 1));
%xlabel('Tempo (s)');
%ylabel('Aceleração (Gs)');
%title('Aceleraçao do Banco');

%line([0.1, 0.2], [35, 35], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
%xline(0.1, 'LineWidth', 2);
%xline(0.2, 'LineWidth', 2);
%ylim([20 75])
%xlim([0.05 0.25])

%legend('Sem controle','Alocação de Polos','LQR','Location','northeast','Orientation','vertical')

%dim = [0.35 0.92 .0001 .0001];
%str1 = 'Zona de Aumento de Lesão';
%annotation('textbox',dim,'String',str1,'FitBoxToText','on');
%dim = [0.55 0.2 .0001 .0001];
%str2 = 'Sem Lesão';
%annotation('textbox',dim,'String',str2,'FitBoxToText','on');

%figure(3)
%plot(t,EntradaDeControle_alocacao, t, EntradaDeControle_lqr)
%legend('Alocação de Polos', 'LQR')
%xlabel('Tempo(s)')
%ylabel('u(t)')
%title('Entrada de controle')

%figure(4)
%plot(t, aceleracoes(:, 1));
%xlabel('Tempo (s)');
%ylabel('Aceleração (Gs)');
%title('Aceleraçao do Banco');

%Controle moderno

%hold on;
% Coordenadas do retângulo (área a ser ampliada)
%x1 = 0.05;
%x2 = 0.25;
%y1 = 25;
%y2 = 75;

% Desenhar o retângulo
%rectangle('Position', [x1, y1, x2-x1, y2-y1], 'EdgeColor', 'r', 'LineWidth', 2);

% Definir as coordenadas do novo eixo (axes) para o gráfico ampliado
%axes('Position', [0.6, 0.6, 0.25, 0.25]);
%box on;

% Desenhar o gráfico ampliado
%plot(t, aceleracoes(:, 1));
%xline(0.1, 'LineWidth', 2);
%xline(0.2, 'LineWidth', 2);
%line([0.1, 0.2], [35, 35], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
%x = [0.5 0.68];
%y = [0.8 0.8];
%annotation('textarrow',x,y,'String','Zona de aumento de lesão')
%x = [0.5 0.7];
%y = [0.3 0.6];
%annotation('textarrow',x,y,'String','Sem lesão')
%xlim([x1, x2]);
%ylim([y1, y2]);

%figure(5)
%plot(t_alocacao, x_alocacao(:,1));
%xlim([0 .8])
%ylim([0.5 2])
%xlabel('Tempo (s)');
%ylabel('Deslocamento (m)');
%title('Deslocamento do Banco por alocação de polos');

%figure(6)
%plot(t_alocacao, aceleracoes_alocacao(:, 1));
%xlabel('Tempo (s)');
%ylabel('Aceleração (Gs)');
%title('Aceleraçao do Banco por alocação de polos');

%figure(7)
%plot(t_lqr, x_lqr(:,1));
%xlim([0 .8])
%ylim([0.5 2])
%xlabel('Tempo (s)');
%ylabel('Deslocamento (m)');
%title('Deslocamento do Banco por LQR');

%figure(8)
%plot(t_alocacao, aceleracoes_lqr(:, 1));
%xlabel('Tempo (s)');
%ylabel('Aceleração (Gs)');
%title('Aceleraçao do Banco por LQR');

%figure('units','centimeters','position',[2 2 18 12],'color','w');
%rlocus(FT1_red);
%title('Lugar das raízes para ganho proporcional Kp^0','fontsize',25);
%xlabel('Eixo Real','fontsize',20);
%ylabel('Eixo Imaginário','fontsize',20);

%set(gca,'FontSize',17);
%grid on;




%Observador


%Alocação de polos

%polosdoobservador

t2 = 0:0.001:2;

e=[0,0,1,1,0,0,0,0];
P_alocacao_obs=pole(sistlqr);
P_al=5*real(P_alocacao_obs)+1i*imag(P_alocacao_obs);




L_aloc_obs = place(A', C1', P_al)';

A_obs = A - L_aloc_obs * C1;
sis_obs = ss(A_obs, B, C1, D1);
% [yobs,t]=impulse(sis_obs);
% yobs = reshape(yobs,size(yobs,1),size(yobs,2),1);



%figure(15);
%pzplot(sis_obs);
%title('Diagrama de Polos do Sistema');
%grid on;
%title('Polos do Observador')
%xlabel('Real');
%ylabel('Imaginário');
%grid on;

[Lambda_obs]=lambda(A, B , C1 ,K_lqr, L_aloc_obs);
% E_obs = [E; zeros(size(A_obs, 1), size(E, 2))];
% sys_obs = ss(Lambda_obs, E_obs, eye(16), zeros(16, 2));
% X01_obs = 0.1*ones(16,1);
% 
% 
% 
%  U_step=ones(length(t),2);
%  
%  % Primeira simulação
% [T1_obs, out1_obs, x1_obs] = lsim(sys_obs, U_step, t, X01_obs);
% 
% figure;
% plot(T1_obs, x1_obs);
% xlabel('Tempo (s)');
% ylabel('Estados do Sistema');
% title('Resposta dos Estados do Sistema a uma Entrada Step');
% grid on;
% legend(arrayfun(@(i) sprintf('x_%d', i), 1:16, 'UniformOutput', false));

% Quarta simulação
%[T1_obs_or, out1_obs_or, x1_obs_or] = lsim(sys_obs, x_lqr, t, X01_obs);

% tss=0.001;
% erro0=zeros(length(t2),16);
% erro0(1,:)=0.1*ones(16,1);
% erro1=vertcat(zeros(8,1), 0.1*ones(8,1))';
% 
% fiii=Lambda_obs;
% phi=expm(0.001*Lambda_obs);
% erroM=erro1;
% 
% for i=2:length(t2)
%     erroM(i,:)=phi*erro1(i-1,:)';
%     erro1(i,:)=[erroM(i,1:8),gradient(erroM(i,1:8))]';
% end


tss=0.001;
erro0=zeros(length(t),16);
erro0(1,:)=0.1*ones(16,1);
erro1=erro0;

fiii=Lambda_obs;
phi=expm(tss*Lambda_obs);

for i=2:length(t)
    erro1(i,:)=phi*erro1(i-1,:)';
end

% 
% figure(123);
%    subplot(3,2,1);
%    plot(t2(1:150),erro1(1:150,2),"linewidth",2);
%    hold on
%    plot(t2(1:150),(erro1(1:150,2)-erro1(1:150,10)),"linewidth",2);
%    grid on
%    title( "Respostas do Deslocamento do Tronco" , FontSize=8);
%    xlabel('Tempo(s)');
%    ylabel('Deslocamento Horizontal(m)');
%    legend('real','estimado');
%    subplot(3,2,2);
%    plot(t2(1:150),erro1(1:150,10),'r',"linewidth",2)
%    grid on
%    title("Erro na Resposta do Deslocamento Tronco", FontSize=8)
%    xlabel('Tempo(s)')
%    ylabel("Deslocamento Horizontal(m)")
%    legend('erro')
%    subplot(3,2,3);
%    plot(t2(1:150) ,erro1(1:150,3),"linewidth",2)
%    hold on
%    plot (t2(1:150) ,erro1(1:150,3)-erro1(1:150,11),"linewidth", 2 )
%    grid on
%    title("Respostas do Deslocamento dos Membros Inferiores", FontSize=8)
%    xlabel("Tempo(s)")
%    ylabel("Deslocamento Horizontal(m)")
%    legend("real", "estimado")
%    subplot(3,2,4);
%    plot(t2(1:150),erro1(1:150,11),"r","linewidth" , 2 )
%    grid on
%    title("Erro na Resposta do Deslocamento dos Membros Inferiores", FontSize=8)
%    xlabel('Tempo(s)')
%    ylabel("Deslocamento Horizontal(m)")
%    legend('erro')
%    subplot(3,2,5);
%    plot(t2(1:150),erro1(1:150,4),"linewidth",2);
%    hold on
%    plot(t2(1:150),(erro1(1:150,4)-erro1(1:150,12)),"linewidth",2);
%    grid on
%    title( "Respostas da Rotação da Cabeça" , FontSize=8);
%    xlabel('Tempo(s)');
%    ylabel('Rotação(Graus)');
%    legend('real','estimado');
%    subplot(3,2,6);
%    plot(t2(1:150),erro1(1:150,12),'r',"linewidth",2)
%    grid on
%    title("Erro na Resposta da Rotação da Cabeça", FontSize=8)
%    xlabel('Tempo(s)')
%    ylabel("Rotação(Graus)")
%    legend('erro')
% 
%    figure(127);
%    subplot(3,2,1);
%    plot(t2(1:150),erro1(1:150,6),"linewidth",2);
%    hold on
%    plot(t2(1:150),(erro1(1:150,6)-erro1(1:150,14)),"linewidth",2);
%    grid on
%    title( "Respostas da Aceleração do Tronco" , FontSize=8);
%    xlabel('Tempo(s)');
%    ylabel('Aceleração Horizontal(m/s^2)');
%    legend('real','estimado');
%    subplot(3,2,2);
%    plot(t2(1:150),erro1(1:150,14),'r',"linewidth",2)
%    grid on
%    title("Erro na Aceleração do Tronco", FontSize=8)
%    xlabel('Tempo(s)')
%    ylabel('Aceleração Horizontal(m/s^2)')
%    legend('erro')
%    subplot(3,2,3);
%    plot(t2(1:150) ,erro1(1:150,7),"linewidth",2)
%    hold on
%    plot (t2(1:150) ,erro1(1:150,7)-erro1(1:150,15),"linewidth", 2 )
%    grid on
%    title("Respostas da Aceleração dos Membros Inferiores", FontSize=8)
%    xlabel("Tempo(s)")
%    ylabel("Aceleração Horizontal(m/s^2)")
%    legend("real", "estimado")
%    subplot(3,2,4);
%    plot(t2(1:150),erro1(1:150,15),"r","linewidth" , 2 )
%    grid on
%    title("Erro na Resposta da Aceleração dos Membros Inferiores", FontSize=8)
%    xlabel('Tempo(s)')
%    ylabel("Aceleração Horizontal(m/s^2)")
%    legend('erro')
%    subplot(3,2,5);
%    plot(t2(1:150),erro1(1:150,8),"linewidth",2);
%    hold on
%    plot(t2(1:150),(erro1(1:150,8)-erro1(1:150,16)),"linewidth",2);
%    grid on
%    title( "Respostas da Aceleração Angular da Cabeça" , FontSize=8);
%    xlabel('Tempo(s)');
%    ylabel('Aceleração Angular(graus/s^2)');
%    legend('real','estimado');
%    subplot(3,2,6);
%    plot(t2(1:150),erro1(1:150,16),'r',"linewidth",2)
%    grid on
%    title("Erro na resposta do Banco", FontSize=8)
%    xlabel('Tempo(s)')
%    ylabel("Aceleração Ângular(graus/s^2)")
%    legend('erro')
% 
% 
% % figure(124);
% %    subplot(2,2,1);
% %    plot(t,erro1(:,1),"linewidth",2);
% %    hold on
% %    plot(t,(erro1(:,1)-erro1(:,9)),"linewidth",2);
% %    grid on
% %    title( "Resposta do Deslocamento do Banco" , FontSize=8);
% %    xlabel('Tempo(s)');
% %    ylabel('Deslocamento Horizontal(m)');
% %    legend('real','estimado');
% %    subplot(2,2,2);
% %    plot(t,erro1(:,9),'r',"linewidth",2)
% %    grid on
% %    title("Erro na resposta do Banco", FontSize=8)
% %    xlabel('Tempo(s)')
% %    ylabel("Deslocamento Horizontal(m)")
% %    legend('erro')
% %    subplot(2,2,3);
% %    plot(t ,erro1(:,5),"linewidth",2)
% %    hold on
% %    plot (t ,erro1(:,5)-erro1(:,13),"linewidth", 2 )
% %    grid on
% %    title("Resposta da Velocidade do Banco", FontSize=8)
% %    xlabel("Tempo(s)")
% %    ylabel("Velocidade (m/s)")
% %    legend("real", "estimado")
% %    subplot(2,2,4);
% %    plot(t,erro1(:,13),"r","linewidth" , 2 )
% %    grid on
% %    title("Erro na Resposta da Velocidade do Banco", FontSize=8)
% %    xlabel('Tempo(s)')
% %    ylabel("Velocidade (m/s)")
% %    legend('erro')
% 
% 
% 
% submatriz = erro1(1:length(t2), [10:12,14:16]);
% 
% 
% 
% %submatriz = erro1(1:length(t), 9:16);
% 
% 
% figure(321);
%     plot(t2(1:150),submatriz(1:150,1:2),'linewidth',2);
%     grid on
%     title("Evolucao do Erro")
%     xlabel("Tempo(s)")
%     ylabel("Deslocamento(m)")
%     legend('$x_t$','$x_p$', 'Interpreter', 'latex')
% figure(322);
%     plot(t2(1:150),submatriz(1:150,3),'linewidth',2);
%     grid on
%     title("Evolucao do Erro")
%     xlabel("Tempo(s)")
%     ylabel("Rotação(°)")
%     legend('$\theta$', 'Interpreter', 'latex')
% figure(323);
%     plot(t2(1:150),submatriz(1:150,4:5),'linewidth',2);
%     grid on
%     title("Evolucao do Erro")
%     xlabel("Tempo(s)")
%     ylabel("Aceleração")
%     legend('$\ddot{x}_t$','$\ddot{x}_p$', 'Interpreter', 'latex')
% figure(324);
%     plot(t2(1:150),submatriz(1:150,6),'linewidth',2);
%     grid on
%     title("Evolucao do Erro")
%     xlabel("Tempo(s)")
%     ylabel("Aceleração Angular(°/s^2)")
%     legend('$\ddot{\theta}$', 'Interpreter', 'latex')

% OBSERVADOR LQ
% 
% 
% Qo = [
%      10000,0,0,0,0,0,0,0;
%      0,10,0,0,0,0,0,0;
%      0,0,1000,0,0,0,0,0;
%      0,0,0,1000,0,0,0,0;
%      0,0,0,0,2500,0,0,0;
%      0,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,1000,0;
%      0,0,0,0,0,0,0,0
% ]; 
% 
% Ro = 0.1; 
% 
% 
% K_obslq=lqr(A,B,Qo,Ro);
% A_obslqr=A-B*K_obslq;
% 
% sis_obslqr=ss(A_obslqr,B,C,D);
% [t_lqrobs,x_lqrobs] = ode45(@(t,x) LQR(t,x,A,B,K_obslq,E), tempo, x0);
% %K_obslq=K_obslq';
% [lambdalq]=lambda(A ,B ,C ,K_lqr ,K_obslq);
% polosobslqr=pole(sis_obslqr);



% erro2=erro0;
% 
% ffi=lambdalq*tss;
% philq1=expm(ffi);
% tss2=0.001;
% philq=expm(lambdalq*tss2);



% for i = 2:length(t)
%     erro2(i,:)=philq*erro2(i-1,:)';
% end

%figure(43);
%pzplot(sis_obslqr);
%title("Polos do Sistema");
%xlabel('Real');
%ylabel('Imaginario');
%grid on
%legend('Polos');

% linhas = min(length(t), 2000);
% submatrizlq=erro2(1:linhas, 9:16);

%figure(188);
%    subplot(2,2,1);
%    plot(t,erro2(:,1),"linewidth",2);
%    hold on
%    plot(t,(erro2(:,1)-erro2(:,9)),"linewidth",2);
%    grid on
%    title( "Resposta do Deslocamento do Banco" , FontSize=8);
%    xlabel('Tempo(s)');
%    ylabel('Deslocamento Horizontal(m)');
%    legend('real','estimado');
%    subplot(2,2,2);
%    plot(t,erro2(:,9),'r',"linewidth",2)
%    grid on
%    title("Erro na resposta do Banco", FontSize=8)
%    xlabel('Tempo(s)')
%    ylabel("Deslocamento Horizontal(m)")
%    legend('erro')
%    subplot(2,2,3);
%    plot(t ,erro2(:,5),"linewidth",2)
%    hold on
%    plot (t ,erro2(:,5)-erro2(:,13),"linewidth", 2 )
%    grid on
%    title("Resposta da Velocidade do Banco", FontSize=8)
%    xlabel("Tempo(s)")
%    ylabel("Velocidade (m/s)")
%    legend("real", "estimado")
%    subplot(2,2,4);
%    plot(t,erro2(:,13),"r","linewidth" , 2 )
%    grid on
%    title("Erro na Resposta da Velocidade do Banco", FontSize=8)
%    xlabel('Tempo(s)')
%    ylabel("Velocidade (m/s)")
%    legend('erro')

    %figure(123);
%    subplot(2,2,1);
%    plot(t,erro1(:,1),"linewidth",2);
%    hold on
%    plot(t,(erro1(:,1)-erro1(:,9)),"linewidth",2);
%    grid on
%    title( "Resposta do Deslocamento do Banco" , FontSize=8);
%    xlabel('Tempo(s)');
%    ylabel('Deslocamento Horizontal(m)');
%    legend('real','estimado');
%    subplot(2,2,2);
%    plot(t,erro1(:,9),'r',"linewidth",2)
%    grid on
%    title("Erro na resposta do Banco", FontSize=8)
%    xlabel('Tempo(s)')
%    ylabel("Deslocamento Horizontal(m)")
%    legend('erro')
%    subplot(2,2,3);
%    plot(t ,erro1(:,5),"linewidth",2)
%    hold on
%    plot (t ,erro1(:,5)-erro1(:,13),"linewidth", 2 )
%    grid on
%    title("Resposta da Velocidade do Banco", FontSize=8)
%    xlabel("Tempo(s)")
%    ylabel("Velocidade (m/s)")
%    legend("real", "estimado")
%    subplot(2,2,4);
%    plot(t,erro1(:,13),"r","linewidth" , 2 )
%    grid on
%    title("Erro na Resposta da Velocidade do Banco", FontSize=8)
%    xlabel('Tempo(s)')
%    ylabel("Velocidade (m/s)")
%    legend('erro')


%figure(37);
%plot(t(1:length(t)),submatrizlq,"linewidth", 2)
%grid on
%title('Evolução do Erro')
%xlabel('Tempo(s)')
%legend('$x_b$','$x_t$','$x_p$','$\theta$','$\dot{x}_b$','$\dot{x}_t$','$\dot{x}_p$','$\dot{\theta}$', 'Interpreter', 'latex')


%Seguidor

h=0.38;
ts=0.01;

[K_seg,sis_seg,yseg_h,t_seg] = seguidordereferenciaconstante(A,B,C4,D4,K_lqr,h,ts);

figure(645784)
    h1=plot(t_seg, yseg_h);
    grid on
    title('Resposta com Seguidor de Referência Constante');
    xlabel('Tempo(s)');
    ylabel('Deslocamento(m)');
    legend('x_b');
    hold on
    h2 = plot(t_seg,0.38*ones(1,length(t_seg)));
    legend([h1,h2],{'Deslocamento do Banco','Referência'})


% Controle Clássico

%FT = tf(sistema(5,1)); 
%G = FT;
% Ganho proporcional ótimo
%Kp = 0.45;

% Função de transferência em malha fechada com Kp
%FT_MF = feedback(Kp * G, 1);

% Função de transferência integrativa
%s = tf('s');
%FT_I = FT_MF / s;

%figure(14)
%rlocus(FT_I);
%grid on 
%ylim([-100 100])
%xlim([-100 0])
%title('Lugar das Raízes para Variação de Ki');

% Ganho integral ótimo
%Ki = 20;

% Função de transferência em malha fechada com Kp e Ki
%FT_MF_comKi = feedback(Kp * G + Ki / s * G, 1);

% Função de transferência integrativa
%s = tf('s');
% Cálculo da função de transferência derivativa
%Kd = 0; % Ganho derivativo
%FT_D = Kd * s;

% Adição da parte derivativa ao sistema em malha fechada
%FT_MF_comKd = feedback(Kp * G + Ki / s * G + Kd * s * G, 1);

%Rootlocusgrafico
%figure(15)
%rlocus(FT_D);
%grid on 
%title('Lugar das Raízes para Variação de Kd');


% Controlador PID
%PI = pid(Kp, Ki, Kd);

% Função de transferência em malha fechada com o controlador PI
%T = feedback(PI * G, 1);
        
%caracteristicas da simulação
%tempoDeSimulacao = 2; %s
%amostra = 0.0001;
%tempo = 0:amostra:tempoDeSimulacao;
%x0 = [0; 0; 0; 0; 0; 0; 0; 0];
%u = 55438.083333333*ones(size(tempo));
%[y, t] = lsim(FT1_red, u, tempo, x0);
%[y_pi, t_pi] = lsim(T, u, tempo, x0);
%y = y/9.81;
%y_pi = y_pi/9.81;

%GraficoDeResposta de maplha fechada pi
%figure(16)
%plot(t, y, t_pi, y_pi)
%line([0.1, 0.2], [35, 35], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
%line([0.1, 0.1], [20, 75], 'Color', 'k', 'LineWidth', 1);
%line([0.2, 0.2], [20, 75], 'Color', 'k', 'LineWidth', 1);
%ylim([20 75])
%xlim([0.05 0.25])
%title('Resposta em malha fechada com controlador PI');
%legend('Malha aberta', 'PI')
%ylabel('Aceleração (Gs)')
%xlabel('Tempo(s)')
%dim = [0.35 0.92 .0001 .0001];
%str1 = 'Zona de Aumento de Lesão';
%annotation('textbox',dim,'String',str1,'FitBoxToText','on');
%dim = [0.55 0.2 .0001 .0001];
%str2 = 'Sem Lesão';
%annotation('textbox',dim,'String',str2,'FitBoxToText','on')


% Ganho proporcional ótimo
%Kp_red = 0.1;

% Função de transferência em malha fechada com Kp
%FT_MF_red = feedback(Kp_red * FT1_red, 1);

% Função de transferência integrativa
%s = tf('s');
%FT_I_red = FT_MF_red / s;

%rootlocuski
%figure(11)
%rlocus(FT_I_red);
%grid on 
%ylim([-100 100])
%xlim([-100 0])
%title('Lugar das Raízes para Variação de Ki');

% Ganho integral ótimo
%Ki_red = 10;

% Função de transferência em malha fechada com Kp e Ki
%FT_MF_comKi_red = feedback(Kp_red * FT1_red + Ki_red / s * FT1_red, 1);

% Função de transferência integrativa
%s = tf('s');
% Cálculo da função de transferência derivativa
%Kd_red = 0; % Ganho derivativo
%FT_D_red = Kd * s;

% Adição da parte derivativa ao sistema em malha fechada
%FT_MF_comKd_red = feedback(Kp_red * FT1_red + Ki_red / s * FT1_red + Kd_red * s * FT1_red, 1);

%rootlocuskd
%figure(12)
%rlocus(FT_D_red);
%grid on 
%title('Lugar das Raízes para Variação de Kd');


% Controlador PID
%PI_red = pid(Kp_red, Ki_red, Kd_red);

% Função de transferência em malha fechada com o controlador PI
%T_red = feedback(PI_red * FT1_red, 1);

%simulacao2
%tempoDeSimulacao = 2; %s
%amostra = 0.0001;
%tempo = 0:amostra:tempoDeSimulacao;
%x0 = [0; 0; 0; 0; 0; 0; 0; 0];
%u = 55438.083333333*ones(size(tempo));
%[y, t] = lsim(FT1_red, u, tempo, x0);
%[y_pi_red, t_pi_red] = lsim(T_red, u, tempo, x0);
%y = y/9.81;
%y_pi_red = y_pi_red/9.81;

%frafico de resposta em malha fechada com controlador pi
%figure(13)
%plot(t, y, t_pi_red, y_pi_red)
%line([0.1, 0.2], [35, 35], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
%line([0.1, 0.1], [20, 75], 'Color', 'k', 'LineWidth', 1);
%line([0.2, 0.2], [20, 75], 'Color', 'k', 'LineWidth', 1);
%ylim([20 75])
%xlim([0.05 0.25])
%title('Resposta em malha fechada com controlador PI');
%legend('Malha aberta', 'PI')
%ylabel('Aceleração (Gs)')
%xlabel('Tempo(s)')
%dim = [0.35 0.92 .0001 .0001];
%str1 = 'Zona de Aumento de Lesão';
%annotation('textbox',dim,'String',str1,'FitBoxToText','on');
%dim = [0.55 0.2 .0001 .0001];
%str2 = 'Sem Lesão';
%annotation('textbox',dim,'String',str2,'FitBoxToText','on')

%grafico comparando pi e pi em ordem reduzida
%figure(17)
%plot(t_pi, y_pi, t_pi_red, y_pi_red)
%ylim([-150 150])
%xlim([0 1])
%title('Comparação PI e PI ordem reduzida');
%legend('PI Original', 'PI - ordem reduzida')
%ylabel('Aceleração (Gs)')
%xlabel('Tempo(s)')
%annotation('textbox',dim,'String',str2,'FitBoxToText','on')

%controlabilidade e observabilidade
%MC = ctrb(A,B);
%n = rank(MC, 0.00001);
%if n == length(A)
%    disp('Controlável')
%else
%    disp('Não controlável')
%end


% MO = obsv(A,C2);
% n = rank(MO, 0.00001);
% if n == length(A)
%     disp('Observável')
% else
%     disp('Não observável')
% end



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

    L = 1.6;
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

    L = 1.6;
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
function dxdt = LQR_Observador(t, x, A, B, C2, K_lqr, L_aloc_obs, E)
    % Separar as variáveis de estado do sistema e do observador
    %n = size(A, 1); % Dimensão do sistema
    % Entrada de referência
    ref = [0; 0; 0; 0; 0; 0; 0; 0]';
    % Função para velocidade na colisão (exemplo)
    r2 = velocidadeNaColisao(t);
    L = 0.38;
    k = 20;
    f = 0.01;
    r1 = @(t) L ./ (1 + exp(-k * (t - f)));
    %r1 = @(t) (0.4339+0.0076+0.0011+0.0002)/ (1 + exp(-50 * (t - 0.036)))-0.0539-0.0076-0.0011-0.0002;
    w = [r1(t); 
        r2
        ];
    
    dx_obs = A*x-L_aloc_obs*C2*x+E*w;
    
    % Concatenar as derivadas
    dxdt = [dx_sys; dx_obs];
end
function[Lambda]=lambda(A,B,C,K,Ko)
Lambda = [A-B*K, B*K;
zeros(size(A)),A-Ko*C];
end
function [K,sys,y,t] = seguidordereferenciaconstante(A,B,C,D,K_ref,h,ts)
    A1=[A  B ; C  D];
    A2=inv(A1);
    No=[zeros(8,1);ones(1,1)];
    Nxu=(A2)*No;
    Nx=Nxu(1:8);
    Nu=Nxu(9);
    K=K_ref*Nx+Nu;
    Be=B*K;
    sys=ss(A-B*K_ref,Be,C,zeros(1,1));
    t=0:ts:2;
    u=stepDataOptions('StepAmplitude',h);
    y=step(sys,t,u);
end
