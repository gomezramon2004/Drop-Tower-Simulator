clear; clc; close;

% Variables
fg = 9.81;                                  % Aceleración por gravedad
m = 0.01;                                   % Masa del dipolo          
mu = 1e6;                                   % Permeabilidad del Metglas 
a = 0.08;                                   % Radio del aro de corriente                  
R = 9e-5;                                   % Resistividad del cobre
p0 = 20;                                    % posicion inicial
v0 = 0;                                     % velocidad inicial

% Constantes
mu_0 = 4*pi*1e-7;                           % Permeabilidad al vacío 
k = (9 * (mu * mu_0).^2 * a.^4) / (4 * R);

% RK4
x0 = 0;                                     % tiempo inicial (x inicial)
y0 = [v0; p0];                              % y inicial (incluye posicon y velocidad)
xf = 2.2;                                     % x final
h = 0.2;                                    % Numero de pasos

% Resuelve
[t_rk4, y_rk4] = rk4(@(t, y) edo(t, y, fg, k, m, a), x0, y0, xf, h);

% Extraer los resultados
v_rk4 = y_rk4(1, :);
p_rk4 = y_rk4(2, :);
aaaaa = diff(v_rk4);
aaaaa(12) = 0;

%----------------------------------Plots----------------------------------%

figure(1)
subplot(1,3,1)
plot(t_rk4, v_rk4, 'b', 'LineWidth', 2);
grid on
xlabel('Tiempo (s)');
ylabel('Velocidad (m/s)');
title('Velocidad en función del tiempo - RK4');

subplot(1,3,2)
plot(t_rk4, p_rk4, 'r', 'LineWidth', 2);
grid on
xlabel('Tiempo (s)');
ylabel('Posicion (C)');
title('Posicion en función del tiempo - RK4');

subplot(1,3,3)
plot(t_rk4, aaaaa, 'r', 'LineWidth', 2);
grid on
xlabel('Tiempo (s)');
ylabel('Aceleración (m/s^2)');
title('Aceleración en función del tiempo - RK4');

%--------------------------------Funciones--------------------------------%

function dydt = edo(t, y, fg, k, m, a)
    % Variables de estado
    v = y(1);
    p = y(2);

    % Derivadas
    vel = -fg - (k / m) * (p.^2 / ((p.^2 + a.^2).^(5/2))) * v;
    pos = v;
    dydt = [vel; pos];
end

function [x, y] = rk4(edo, x0, y0, xf, h)
    x = x0;
    y = y0;
    i = 1;
    while i * h <= xf
        k1 = edo(x(i), y(:, i));
        k2 = edo(x(i) + 0.5*h, y(:, i) + 0.5*k1*h);
        k3 = edo(x(i) + 0.5*h, y(:, i) + 0.5*k2*h);
        k4 = edo(x(i) + h, y(:, i) + k3*h);
        x(i+1) = x(i) + h;
        y(:, i+1) = y(:, i) + (k1 + 2*k2 + 2*k3 + k4)*h/6;
        i=i+1;
    end
end
