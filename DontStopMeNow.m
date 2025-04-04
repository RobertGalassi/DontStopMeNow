%% pulizia iniziale e impostazioni di default
close all; clc; clear;
set(0, 'defaultFigureColor', 'k');   % finestra con sfondo nero
set(0, 'defaultAxesColor', 'k');       % assi con sfondo nero

%% parametri fisici
G = 6.67430e-11;             % [m^3/(kg*s^2)] costante di gravitá
M = 1.989e30;                % [kg] massa del sole
mTerra = 5.972e24;           % [kg] massa della terra
L = 2.53e37;                 % [kg*m^2/s] momento angolare

%% simulazione dell'orbita instabile (spirale)
tempoInizio = 0;
tempoFine = 100000 * 365.25 * 24 * 3600;  % 100.000 anni in secondi
numfotoinstabile = 1000;                  % numero di fotogrammi instabile
tempoVec = linspace(tempoInizio, tempoFine, numfotoinstabile);
tempoAnni = tempoVec / (3600 * 24 * 365.25);  % tempo in anni

r0Instabile = 1.826e10;  % [m] condizione iniziale per l'orbita instabile
funDiff = @(t, r) sqrt((4 * G * M * mTerra) / (3 * r^(3/2)) - (L^2)/(mTerra^2 * r^2));
opts = odeset('RelTol',1e-8, 'AbsTol',1e-10);
[tempoSol, rSol] = ode45(funDiff, tempoVec, r0Instabile, opts);
if any(rSol <= 0)
    error('il raggio è negativo o nullo. verifica le equazioni.');
end

%% calcolo dell'orbita chiusa (circolare)
% r_chiusa = L^4/(G^2 * M^2 * mTerra^4)
rChiusa = L^4 / (G^2 * M^2 * mTerra^4);
fprintf('raggio orbita chiusa: %e m\n', rChiusa);
vChiusa = L / (mTerra * rChiusa);   % v = L/(mTerra*rChiusa)
T = 2*pi * mTerra * rChiusa^2 / L;     % periodo orbitale

numfotoorbchiusa = 500;               % numero di fotogrammi per orbita chiusa
tempiOrbchiusa = linspace(0, T, numfotoorbchiusa);
phiOrbchiusa = 2*pi * (tempiOrbchiusa/T);
xOrbchiusa = rChiusa * cos(phiOrbchiusa);
yOrbchiusa = rChiusa * sin(phiOrbchiusa);

%% calcolo dell'energia meccanica nell'orbita chiusa
term1 = L^2 / (2 * mTerra * rChiusa^2);
term2 = 2 * G * M * mTerra / (3 * rChiusa^(3/2));
energiaMechOrbchiusa = term1 - term2;

%% impostazione della figura con 3 subplot
figure('Position', [100, 100, 1600, 900]);

%% subplot 1: grafico dell'evoluzione del raggio (instabile)
asse1 = subplot(1,3,1);
set(asse1, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w');
grid(asse1, 'on'); hold(asse1, 'on');
plot(asse1, tempoAnni, rSol, 'w-', 'LineWidth', 1.5);
xlabel(asse1, 'tempo [anni]', 'Color', 'w');
ylabel(asse1, 'raggio [m]', 'Color', 'w');
title(asse1, 'evoluzione raggio instabile', 'Color', 'w');
rmin = min(rSol); rmax = max(rSol);
lineaTempo = line(asse1, [tempoAnni(1) tempoAnni(1)], [rmin, rmax], 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
markerTempo = plot(asse1, tempoAnni(1), rSol(1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

%% subplot 2: animazione dell'orbita instabile (spirale)
asse2 = subplot(1,3,2);
set(asse2, 'Color','k','XColor','w','YColor','w'); axis(asse2, 'equal'); hold(asse2, 'on');
rmaxPlotInst = max(rSol)*1.1;
axis(asse2, [-rmaxPlotInst, rmaxPlotInst, -rmaxPlotInst, rmaxPlotInst]);
title(asse2, 'animazione orbita instabile', 'Color', 'w');
plot(asse2, 0, 0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'y');  % stella centrale

rPianetaInst = rmaxPlotInst * 0.02;  % dimensione del pianeta instabile
angPianeta = linspace(0, 2*pi, 50);
faseInst = 0;      % fase iniziale
dFaseInst = 0.1;   % incremento di fase per l'animazione
xPianetaI0 = rSol(1) * cos(faseInst);
yPianetaI0 = rSol(1) * sin(faseInst);
xPianetaInst = xPianetaI0 + rPianetaInst * cos(angPianeta);
yPianetaInst = yPianetaI0 + rPianetaInst * sin(angPianeta);
patchInstabile = patch(asse2, xPianetaInst, yPianetaInst, 'b', 'EdgeColor', 'b');
tracciaXInst = [];
tracciaYInst = [];
tracciaInstabile = plot(asse2, NaN, NaN, 'w-', 'LineWidth', 1);

%% subplot 3: animazione dell'orbita chiusa (stabile)
asse3 = subplot(1,3,3);
set(asse3, 'Color', 'k','XColor','w','YColor','w'); axis(asse3, 'equal'); hold(asse3, 'on');
rmaxPlotChiusa = rChiusa * 1.2;
axis(asse3, [-rmaxPlotChiusa, rmaxPlotChiusa, -rmaxPlotChiusa, rmaxPlotChiusa]);
title(asse3, 'animazione orbita chiusa', 'Color', 'w');
plot(asse3, 0, 0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'y');  % stella centrale
angCerchio = linspace(0, 2*pi, 200);
plot(asse3, rChiusa * cos(angCerchio), rChiusa * sin(angCerchio), 'w--', 'LineWidth', 1.5);
rPianetaChiusa = rmaxPlotChiusa * 0.02;
% definisco il colore azzurro
coloreAzzurro = [0, 0.749, 1];
xPianetaC0 = xOrbchiusa(1) + rPianetaChiusa * cos(angPianeta);
yPianetaC0 = yOrbchiusa(1) + rPianetaChiusa * sin(angPianeta);
patchChiusa = patch(asse3, xPianetaC0, yPianetaC0, coloreAzzurro, 'EdgeColor', coloreAzzurro);
text(asse3, -rmaxPlotChiusa*0.9, rmaxPlotChiusa*0.9, sprintf('energia meccanica = %e J', energiaMechOrbchiusa), 'Color', 'w', 'FontSize', 10, 'BackgroundColor', 'k');

%% loop di animazione sincronizzato
numFotoInst = length(tempoSol);
NfotoTotali = max(numFotoInst, numfotoorbchiusa);
for fot = 1:NfotoTotali
    %% aggiornamento asse1: linea e marker sul grafico del tempo
    idxTempo = min(fot, numFotoInst);
    tempoCorr = tempoAnni(idxTempo);
    set(lineaTempo, 'XData', [tempoCorr, tempoCorr]);
    set(markerTempo, 'XData', tempoCorr, 'YData', rSol(idxTempo));
    
    %% aggiornamento asse2: animazione orbita instabile
    idxInst = min(fot, numFotoInst);
    rAttInst = rSol(idxInst);
    % aggiorno la fase per far ruotare il pianeta instabile
    faseInst = faseInst + dFaseInst;
    xAttInst = rAttInst * cos(faseInst);
    yAttInst = rAttInst * sin(faseInst);
    xPianetaInst = xAttInst + rPianetaInst * cos(angPianeta);
    yPianetaInst = yAttInst + rPianetaInst * sin(angPianeta);
    set(patchInstabile, 'XData', xPianetaInst, 'YData', yPianetaInst);
    % aggiorno la traccia
    tracciaXInst(end+1) = xAttInst;
    tracciaYInst(end+1) = yAttInst;
    set(tracciaInstabile, 'XData', tracciaXInst, 'YData', tracciaYInst);
    
    %% aggiornamento asse3: animazione orbita chiusa
    idxChiusa = round((fot-1) / (NfotoTotali-1) * (numfotoorbchiusa-1)) + 1;
    xAttChiusa = xOrbchiusa(idxChiusa);
    yAttChiusa = yOrbchiusa(idxChiusa);
    xPianetaChiusa = xAttChiusa + rPianetaChiusa * cos(angPianeta);
    yPianetaChiusa = yAttChiusa + rPianetaChiusa * sin(angPianeta);
    set(patchChiusa, 'XData', xPianetaChiusa, 'YData', yPianetaChiusa);
    
    drawnow;
    pause(0.02);  % breve pausa per fluidità
end

disp('animazione completata.');
