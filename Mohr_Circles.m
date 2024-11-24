clear
clc

sigma_x = input("Sigma x: (MPa) >> "); % Sforzo normale nella direzione x
sigma_y = input("Sigma y: (MPa) >> ") ; % Sforzo normale nella direzione y
sigma_z = input("Sigma z: (MPa) >> ") ; % Sforzo normale nella direzione z
tau_xy = input("Tau xy: (MPa) >> ");  % Sollecitazione tangenziale (di taglio) nella direzione xy
tau_xz = input("Tau zx: (MPa) >> ");    % Sollecitazione tangenziale (di taglio) nella direzione xz
tau_yz = input("Tau yz: (MPa) >> ");    % Sollecitazione tangenziale (di taglio) nella direzione yz

clc

% Inizializzo il piano dello sforzo piano
plane = '';
theta_deg = 0;

% Creazione della matrice del tensore degli sforzi
sigma_tensor = [sigma_x, tau_xy, tau_xz;
       tau_xy, sigma_y, tau_yz;
       tau_xz, tau_yz, sigma_z];

if ~(sigma_tensor(2) == 0)
    % lo stato di sforzo piano è nel piano xy
    plane = 'xy';
else 
    if ~(sigma_tensor(3) == 0)
    % lo stato di sforzo piano è nel piano xz
    plane = 'xz';
    else
    % lo stato di sforzo piano è nel piano yz
    plane = 'yz';
    end
end

% Calcolare gli autovalori e gli autovettori
[autovettori, autovalori] = eig(sigma_tensor);

% Gli autovalori rappresentano gli sforzi principali
[sforzi_principali, indices] = sort(diag(autovalori), 'descend');

% Gli autovettori rappresentano le direzioni principali
direzioni_principali = autovettori(indices,:);


sigma_1 = sforzi_principali(1);
sigma_2 = sforzi_principali(2);
sigma_3 = sforzi_principali(3);
n_1 = direzioni_principali(1,:);
n_2 = direzioni_principali(2,:);
n_3 = direzioni_principali(3,:);
tau_2_3 = sort([-0.5*(sigma_2 - sigma_3),0.5*(sigma_2 - sigma_3)], 'descend');
tau_1_3 = sort([-0.5*(sigma_1 - sigma_3),0.5*(sigma_1 - sigma_3)], 'descend');
tau_1_2 = sort([-0.5*(sigma_1 - sigma_2),0.5*(sigma_1 - sigma_2)], 'descend');

% Calcolare il centro e il raggio del cerchio di Mohr
center_1 = (sigma_1 + sigma_2)/ 2;
radius_1 = (sigma_1 - sigma_2)/2;
center_2 = (sigma_2 + sigma_3)/ 2;
radius_2 = (sigma_3 - sigma_2)/2;
center_3 = (sigma_1 + sigma_3)/ 2;
radius_3 = (sigma_1 - sigma_3)/2;

% Creare i punti per il cerchio di Mohr
theta = linspace(0, 2*pi, 100);  % Angoli da 0 a 2*pi
%sigma = center + radius * cos(2*theta);  % Componenti della tensione normale
%tau = radius * sin(2*theta);  % Componenti della tensione di taglio
sigma = (sigma_x + sigma_y)/2 +((sigma_x - sigma_y)/2)*cos(2*theta) + tau_xy*sin(2*theta);
tau = (-(sigma_x - sigma_y)/2)*sin(2*theta) + tau_xy*cos(2*theta);

% Fattore di proporsionalità assiale
fpa = max([radius_1, radius_2, radius_3])/3;
r = fpa;

% Display delle variabili
fprintf('----------------Mohr Circle v1.0 by AV----------------\n\n Stress tensor:\n');
disp(sigma_tensor);
fprintf('Stress I: %.2f\n',sigma_1);
fprintf('Stress II: %.2f\n',sigma_2);
fprintf('Stress III: %.2f\n\n',sigma_3);
fprintf('N I: [%.3f, %.3f, %.3f]\n',n_1(1), n_1(2),n_1(3));
fprintf('N II: [%.3f, %.3f, %.3f]\n',n_2(1), n_2(2),n_2(3));
fprintf('N III: [%.3f, %.3f, %.3f]\n',n_3(1), n_3(2),n_3(3));


% Disegnare il cerchio di Mohr
if nnz(sigma_tensor == 0) - nnz([sigma_x, sigma_y, sigma_z] == 0) >= 4
    % lo stato di sforzo è piano

figure;
plot(center_1 + radius_1 * cos(2*theta),radius_1 * sin(2*theta) , 'm-', 'LineWidth', 1);
hold on
plot(center_2 + radius_2 * cos(2*theta),radius_2 * sin(2*theta) , 'm-', 'LineWidth', 1);
plot(center_3 + radius_3 * cos(2*theta),radius_3 * sin(2*theta) , 'b-', 'LineWidth', 1);
xlabel('Sforzo normale (MPa)');
ylabel('Sforzo tangenziale (MPa)');

% Tracciare i punti di sforzo su Mohr

text(sigma_1, 0, sprintf('\\sigma_1 = %.2f MPa', sigma_1), 'VerticalAlignment', 'bottom');
text(sigma_2, 0, sprintf('\\sigma_2 = %.2f MPa', sigma_2), 'VerticalAlignment', 'top');
text(sigma_3, 0, sprintf('\\sigma_3 = %.2f MPa', sigma_3), 'VerticalAlignment', 'top');

if isequal(plane,'xy')
theta_deg = 0.5*rad2deg(atan((2*tau_xy)/(sigma_x - sigma_y)));
theta_rad = deg2rad(theta_deg);

% Calcolare le coordinate dei punti sull'arco
theta_vals = linspace(-2*theta_rad, 0, 100);   % Creiamo 100 punti tra 0 e l'angolo theta
x = r * cos(theta_vals) + (sigma_y + sigma_x)/2;    % Coordinate x dei punti sull'arco
y = r * sin(theta_vals);
plot(x, y, '--r', 'LineWidth', 1);
text(sigma_y, tau_xy, sprintf('A (\\sigma_y %.2f MPa, \\tau_xy %.2f MPa)', sigma_y, tau_xy), 'VerticalAlignment', 'top');
text(sigma_x, -tau_xy, sprintf('B (\\sigma_x %.2f MPa, \\tau_xy %.2f MPa)', sigma_x, tau_xy), 'VerticalAlignment', 'top');

% Disegnare le linee che collegano i punti (sigma_x, -tau_xy) e (sigma_y, tau_xy)
plot([sigma_x, sigma_y], [-tau_xy, tau_xy], 'b--', 'LineWidth', 1);
plot(sigma_x, -tau_xy, 'b*', 'MarkerSize', 5, 'LineWidth', 2);  
plot(sigma_y, tau_xy, 'b*', 'MarkerSize', 5, 'LineWidth', 2);
end

if isequal(plane,'xz')
theta_deg = 0.5*rad2deg(atan((2*tau_xz)/(sigma_x - sigma_z)));
theta_rad = deg2rad(theta_deg);

% Calcolare le coordinate dei punti sull'arco
theta_vals = linspace(-2*theta_rad, 0, 100);   % Creiamo 100 punti tra 0 e l'angolo theta
x = r * cos(theta_vals) + (sigma_z + sigma_x)/2;    % Coordinate x dei punti sull'arco
y = r * sin(theta_vals);
plot(x, y, '--r', 'LineWidth', 1);

text(sigma_x, tau_xz, sprintf('A (\\sigma_x %.2f MPa, \\tau_xy %.2f MPa)', sigma_x, tau_xz), 'VerticalAlignment', 'top');
text(sigma_z, -tau_xz, sprintf('B (\\sigma_z %.2f MPa, \\tau_xy %.2f MPa)', sigma_z, tau_xz), 'VerticalAlignment', 'top');

% Disegnare le linee che collegano i punti 
plot([sigma_x, sigma_z], [-tau_xz, tau_xz], 'b--', 'LineWidth', 1);
plot(sigma_x, -tau_xz, 'b*', 'MarkerSize', 5, 'LineWidth', 2);  
plot(sigma_z, tau_xz, 'b*', 'MarkerSize', 5, 'LineWidth', 2);
end

if isequal(plane,'yz')
theta_deg = 0.5*rad2deg(atan((2*tau_xz)/(sigma_x - sigma_z)));
theta_rad = deg2rad(theta_deg);

% Calcolare le coordinate dei punti sull'arco
theta_vals = linspace(-2*theta_rad, 0, 100);   % Creiamo 100 punti tra 0 e l'angolo theta
x = r * cos(theta_vals) + (sigma_y + sigma_z)/2;    % Coordinate x dei punti sull'arco
y = r * sin(theta_vals);
plot(x, y, '--r', 'LineWidth', 1);
text(sigma_y, -tau_yz, sprintf('A (\\sigma_y %.2f MPa, \\tau_yz %.2f MPa)', sigma_x, tau_yz), 'VerticalAlignment', 'top');
text(sigma_z, tau_yz, sprintf('B (\\sigma_z %.2f MPa, \\tau_yz %.2f MPa)', sigma_z, tau_yz), 'VerticalAlignment', 'top');

% Disegnare le linee che collegano i punti
plot([sigma_y, sigma_z], [-tau_yz, tau_yz], 'b--', 'LineWidth', 1);
plot(sigma_y, -tau_yz, 'b*', 'MarkerSize', 5, 'LineWidth', 2);  
plot(sigma_z, tau_yz, 'b*', 'MarkerSize', 5, 'LineWidth', 2);
end


% Tracciare i punti degli sforzi principali
plot(sigma_1, 0, 'k*', 'MarkerSize', 5, 'LineWidth', 2); % Punto sigma_1
plot(sigma_2, 0, 'k*', 'MarkerSize', 5, 'LineWidth', 2); % Punto sigma_2
plot(sigma_3, 0, 'k*', 'MarkerSize', 5, 'LineWidth', 2); % Punto sigma_3

% Tracciare i tau principali
plot((sigma_2 + sigma_1)/2, tau_1_2(1), 'b*', 'MarkerSize', 5, 'LineWidth', 2); % Punto tau 12
plot((sigma_1 + sigma_3)/2, tau_1_3(1), 'b*', 'MarkerSize', 5, 'LineWidth', 2); % Punto tau 13
plot((sigma_2 + sigma_3)/2, tau_2_3(1), 'b*', 'MarkerSize', 5, 'LineWidth', 2); % Punto tau 23

text((sigma_2 + sigma_1)/2, tau_1_2(1), sprintf('\\tau_1 = %.1f MPa', tau_1_2(1)), 'VerticalAlignment', 'bottom');
text((sigma_1 + sigma_3)/2, tau_1_3(1), sprintf('\\tau_2 = %.1f MPa', tau_1_3(1)), 'VerticalAlignment', 'bottom');
text((sigma_2 + sigma_3)/2, tau_2_3(1), sprintf('\\tau_3 = %.1f MPa', tau_2_3(1)), 'VerticalAlignment', 'bottom');


%Asse x
plot([min([sigma_3,0]) - fpa, max([sigma_1, 0])+ fpa], [0, 0], 'k', 'LineWidth', 1);
%quiver(max([sigma_1, 0])+ fpa, 0, 5, 0, 'LineWidth', 1, 'Color', 'k'); % Freccia sull'asse X
text(max([sigma_1, 0])+ fpa, 0, sprintf('\\sigma'),'FontSize', 12, 'VerticalAlignment', 'bottom');
%Asse y
plot([0, 0],[-max([radius_1, radius_2, radius_3]) - fpa, max([radius_1, radius_2, radius_3]) + fpa], 'k', 'LineWidth', 1);
%quiver(0, max([radius_1, radius_2, radius_3]) + fpa, 0, 5, 'LineWidth', 1, 'Color', 'k'); % Freccia sull'asse Y
text(0,max([radius_1, radius_2, radius_3]) + fpa, sprintf('\\tau'),'FontSize', 10, 'VerticalAlignment', 'bottom');
% Asse X (da -1 a 1)
%plot([-300, 300], [0, 0], 'k', 'LineWidth', 1); % 'k' è il colore nero
% Asse Y (da -1 a 1)
%plot([0,0], [-300, 300], 'k', 'LineWidth', 1); % 'k' è il colore nero

% Titolo e legenda
title('Mohr Circle');


% Impostazioni grafiche
%axis([-200 200 -200 200]); 
axis equal
grid on;
else
    fprintf("\nLo stato di sforzo non è piano\n")
end