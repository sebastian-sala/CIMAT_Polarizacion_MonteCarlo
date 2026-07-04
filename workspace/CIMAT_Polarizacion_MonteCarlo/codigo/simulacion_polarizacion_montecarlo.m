%% Simulacion Monte Carlo de polarizacion atmosferica
% Modelo de transporte radiativo polarizado con formalismo Stokes-Mueller.
% El codigo simula fotones como caminantes aleatorios en una atmosfera
% plano-paralela y calcula la elipticidad neta V/I de los fotones
% transmitidos con dispersion multiple.

clear; clc; close all;

%% Parametros iniciales
N_fotones = 10000;
tau_max = 2.0;
S0 = [1; 0; 0; 0]; % Luz solar inicialmente no polarizada [I; Q; U; V]

Stokes_finales = [];
colisiones_totales = [];

h = waitbar(0, 'Simulando fotones...');

%% Caminata aleatoria
for i = 1:N_fotones
    z = tau_max;   % posicion inicial en profundidad optica
    mu_z = -1;     % direccion inicial hacia el suelo
    S = S0;
    n_colisiones = 0;

    while (z > 0) && (z < tau_max + 1e-6)
        % Vuelo libre: inversion de la distribucion acumulativa
        r1 = rand();
        delta_tau = -log(r1);
        z = z + delta_tau * mu_z;

        if z <= 0
            Stokes_finales = [Stokes_finales, S]; %#ok<AGROW>
            colisiones_totales = [colisiones_totales, n_colisiones]; %#ok<AGROW>
            break;
        elseif z >= tau_max
            break;
        end

        % Evento de dispersion
        n_colisiones = n_colisiones + 1;
        phi = 2*pi*rand();

        aceptado = false;
        while ~aceptado
            theta_candidato = pi*rand();
            P_rayleigh = 0.75*(1 + cos(theta_candidato)^2);
            if rand()*1.5 < P_rayleigh
                theta = theta_candidato;
                aceptado = true;
            end
        end

        mu_z = mu_z*cos(theta) - sqrt(max(0, 1 - mu_z^2))*sin(theta)*cos(phi);
        S = MatrizAtmosferaReal(theta) * MatrizRotacion(phi) * S;
    end

    if mod(i, 1000) == 0
        waitbar(i/N_fotones, h);
    end
end

close(h);

%% Resultados
idx_multiple = colisiones_totales >= 2;
Stokes_multi = Stokes_finales(:, idx_multiple);
S_promedio = mean(Stokes_multi, 2);
Elipticidad = S_promedio(4) / S_promedio(1);

fprintf('Simulacion completada.\\n');
fprintf('Fotones totales detectados en el suelo: %d\\n', size(Stokes_finales, 2));
fprintf('Fotones con dispersion multiple: %d\\n', sum(idx_multiple));
fprintf('Vector de Stokes Promedio [I, Q, U, V]: [%.4e, %.4e, %.4e, %.4e]\\n', ...
    S_promedio(1), S_promedio(2), S_promedio(3), S_promedio(4));
fprintf('Elipticidad Neta (V/I): %.4e\\n', Elipticidad);

%% Visualizacion
figure('Name', 'Resultados de polarizacion atmosferica', ...
    'Position', [100, 100, 1200, 400]);

% Histograma de colisiones
subplot(1, 3, 1);
histogram(colisiones_totales, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'w');
title('Distribucion de Colisiones (Monte Carlo)');
xlabel('Numero de eventos de dispersion');
ylabel('Frecuencia (fotones)');
grid on;
xlim([0 max(colisiones_totales)]);

% Esfera de Poincare
subplot(1, 3, 2);
[x_sph, y_sph, z_sph] = sphere(40);
surf(x_sph, y_sph, z_sph, 'FaceAlpha', 0.1, 'EdgeColor', [0.8 0.8 0.8]);
hold on;
line([-1.2 1.2], [0 0], [0 0], 'Color', 'k', 'LineStyle', '--');
line([0 0], [-1.2 1.2], [0 0], 'Color', 'k', 'LineStyle', '--');
line([0 0], [0 0], [-1.2 1.2], 'Color', 'k', 'LineStyle', '--');

Ip = sqrt(S_promedio(2)^2 + S_promedio(3)^2 + S_promedio(4)^2);
s1 = S_promedio(2) / Ip;
s2 = S_promedio(3) / Ip;
s3 = S_promedio(4) / Ip;

quiver3(0, 0, 0, s1, s2, s3, 0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
scatter3(s1, s2, s3, 100, 'r', 'filled', 'MarkerEdgeColor', 'k');
title('Estado en la Esfera de Poincare');
xlabel('Q (lineal H/V)');
ylabel('U (lineal +/-45 deg)');
zlabel('V (circular)');
axis equal; view(135, 30); grid on; hold off;

% Elipse de polarizacion
subplot(1, 3, 3);
psi = 0.5 * atan2(S_promedio(3), S_promedio(2));
chi = 0.5 * asin(S_promedio(4) / Ip);

a = sqrt(Ip) * cos(chi);
b = sqrt(Ip) * sin(chi);

t = linspace(0, 2*pi, 100);
E_x_prime = a * cos(t);
E_y_prime = b * sin(t);

E_x = E_x_prime*cos(psi) - E_y_prime*sin(psi);
E_y = E_x_prime*sin(psi) + E_y_prime*cos(psi);

plot(E_x, E_y, 'r', 'LineWidth', 2);
hold on;
line([-1.1*a 1.1*a], [0 0], 'Color', 'k', 'LineStyle', ':');
line([0 0], [-1.1*a 1.1*a], 'Color', 'k', 'LineStyle', ':');
plot([-a*cos(psi) a*cos(psi)], [-a*sin(psi) a*sin(psi)], 'b--');
title(sprintf('Elipse de Polarizacion (V/I = %.1e)', Elipticidad));
xlabel('E_x');
ylabel('E_y');
axis equal; grid on;
xlim([-1.2*a 1.2*a]); ylim([-1.2*a 1.2*a]);
hold off;

%% Funciones auxiliares
function M = MatrizAtmosferaReal(theta)
    c = cos(theta);
    c2 = c^2;
    s = sin(theta);

    % Matriz de Rayleigh pura
    M_ray = (3/4) * [ ...
        c2 + 1, c2 - 1, 0, 0; ...
        c2 - 1, c2 + 1, 0, 0; ...
        0, 0, 2*c, 0; ...
        0, 0, 0, 2*c ...
    ];

    % Matriz fenomenologica de aerosoles tipo Mie simplificada
    M_mie = [ ...
        1 + c2, c2 - 1, 0, 0; ...
        c2 - 1, 1 + c2, 0, 0; ...
        0, 0, 2*c, 0.5*s; ...
        0, 0, -0.5*s, 2*c ...
    ];

    M = 0.90 * M_ray + 0.10 * M_mie;
end

function L = MatrizRotacion(phi)
    c2p = cos(2*phi);
    s2p = sin(2*phi);

    L = [ ...
        1, 0, 0, 0; ...
        0, c2p, s2p, 0; ...
        0, -s2p, c2p, 0; ...
        0, 0, 0, 1 ...
    ];
end
