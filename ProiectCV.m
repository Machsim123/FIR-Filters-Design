clc
clear all
close all

%% date student

ng = 1;
ns = 1;

%% apelam fct

[omega_p,omega_s,M] = PS_PRJ_2_Faza_1ab(ng,ns);

%% memoram datele pentru mai tarziu, eventual

omega_pp = omega_p;
omega_ss = omega_s;
M_ajustat_ab = M - 1;
Mm = M_ajustat_ab;

%% a

% mat 1: caract spectrale {omega_p, omega_s}
matrix1 = figure('Name', 'Faza 1 a - Matricea 1: Caracteristici spectrale {omega_p, omega_s}');
for i = 1:length(M_ajustat_ab) % stim deja ca e 8 - 1 = 7
    % definirea vect de frecvente si amplitudini pentru firls
    W = [0 omega_p/pi omega_s/pi 1];
    A = [1 1 0 0];
    
    h_firls = firls(M_ajustat_ab(i), W, A); % proiectam filtrul cu firls
    [H, f] = freqz(h_firls, 1, 5000); % calc rasp in frecv

    % linia 1: spectrul in coord liniare
    subplot(3, length(M_ajustat_ab), i);
    plot(f/pi, abs(H));
    hold on;
    xline(omega_p/pi, '--r');
    xline(omega_s/pi, '--r');
    title(sprintf('LS: M = %d', M_ajustat_ab(i)));
    xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
    grid on;

    % linia 2: spectrul in dB
    subplot(3, length(M_ajustat_ab), i + length(M_ajustat_ab));
    plot(f/pi, 20*log10(abs(H)));
    hold on;
    xline(omega_p/pi, '--r');
    xline(omega_s/pi, '--r');
    xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})| [dB]');
    grid on;

    % linia 3: faza
    subplot(3, length(M_ajustat_ab), i + 2*length(M_ajustat_ab));
    plot(f/pi, unwrap(angle(H)));
    xlabel('\omega / \pi'); ylabel('Fază [rad]');
    grid on;
end

% mat 2: caract spectrale {pi - omega_s, pi - omega_p}
matrix2 = figure('Name', 'Faza 1 a - Matricea 2: Caracteristici spectrale {pi - omega_s, pi - omega_p}');
for i = 1:length(M_ajustat_ab)
    % transformare pulsatii
    omega_p_prim = pi - omega_s;
    omega_s_prim = pi - omega_p;
    W = [0 omega_p_prim/pi omega_s_prim/pi 1];

    h_firls = firls(M_ajustat_ab(i), W, A); % proiectam filtrul cu firls

    [H, f] = freqz(h_firls, 1, 5000); % calc rasp in frecv

    % linia 1: spectrul in coord liniare
    subplot(3, length(M_ajustat_ab), i);
    plot(f/pi, abs(H));
    hold on;
    xline(omega_p_prim/pi, '--r');
    xline(omega_s_prim/pi, '--r');
    title(sprintf('LS: M = %d', M_ajustat_ab(i)));
    xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
    grid on;

    % linia 2: spectrul in dB
    subplot(3, length(M_ajustat_ab), i + length(M_ajustat_ab));
    plot(f/pi, 20*log10(abs(H)));
    hold on;
    xline(omega_p_prim/pi, '--r');
    xline(omega_s_prim/pi, '--r');
    xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})| [dB]');
    grid on;

    % linia 3: faza
    subplot(3, length(M_ajustat_ab), i + 2*length(M_ajustat_ab));
    plot(f/pi, unwrap(angle(H)));
    xlabel('\omega / \pi'); ylabel('Fază [rad]');
    grid on;
end

% mat 3: secv pondere
matrix3 = figure('Name', 'Faza 1 a - Matricea 3: Secvențele pondere');
for i = 1:length(M_ajustat_ab)
    % proiectarea pt {omega_p, omega_s}
    W = [0 omega_p/pi omega_s/pi 1];
    h1 = firls(M_ajustat_ab(i), W, A);

    % proiectarea pt {pi - omega_s, pi - omega_p}
    omega_p_prim = pi - omega_s;
    omega_s_prim = pi - omega_p;
    W_prime = [0 omega_p_prim/pi omega_s_prim/pi 1];
    h2 = firls(M_ajustat_ab(i), W_prime, A);

    % graficul secv pondere
    subplot(length(M_ajustat_ab), 2, 2*i - 1);
    stem(h1, 'filled');
    title(sprintf('Secv. pondere LS: M = %d, {omega_p, omega_s}', M_ajustat_ab(i)));
    xlabel('n'); ylabel('h[n]');
    grid on;

    subplot(length(M_ajustat_ab), 2, 2*i);
    stem(h2, 'filled');
    title(sprintf('Secv. pondere LS: M = %d, {pi-omega_s, pi-omega_p}', M_ajustat_ab(i)));
    xlabel('n'); ylabel('h[n]');
    grid on;
end

%% b

% mat 1: caract spectrale {omega_p, omega_s}
matrix1 = figure('Name', 'Faza 1 b - Matricea 1: Caracteristici spectrale {omega_p, omega_s}');
for i = 1:length(M_ajustat_ab)
    % def vect de frecv si amplitudini pt firpm
    W = [0 omega_p/pi omega_s/pi 1];
    A = [1 1 0 0];

    h_firpm = firpm(M_ajustat_ab(i), W, A); % proiectarea filtrului - firpm

    [H, f] = freqz(h_firpm, 1, 5000); % calc rasp in frecv

    % linia 1: Spectrul in coord liniare
    subplot(3, length(M_ajustat_ab), i);
    plot(f/pi, abs(H));
    hold on;
    xline(omega_p/pi, '--r');
    xline(omega_s/pi, '--r');
    title(sprintf('PM: M = %d', M_ajustat_ab(i)));
    xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
    grid on;

    % linia 2: spectrul in dB
    subplot(3, length(M_ajustat_ab), i + length(M_ajustat_ab));
    plot(f/pi, 20*log10(abs(H)));
    hold on;
    xline(omega_p/pi, '--r');
    xline(omega_s/pi, '--r');
    xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})| [dB]');
    grid on;

    % linia 3: faza
    subplot(3, length(M_ajustat_ab), i + 2*length(M_ajustat_ab));
    plot(f/pi, unwrap(angle(H)));
    xlabel('\omega / \pi'); ylabel('Fază [rad]');
    grid on;
end

% mat 2: caract spectrale pt {pi - omega_s, pi - omega_p}
matrix2 = figure('Name', 'Faza 1 b - Matricea 2: Caracteristici spectrale {pi - omega_s, pi - omega_p}');
for i = 1:length(M_ajustat_ab)
    % transf pulsatii
    omega_p_prime = pi - omega_s;
    omega_s_prime = pi - omega_p;
    W = [0 omega_p_prime/pi omega_s_prime/pi 1];

    h_firpm = firpm(M_ajustat_ab(i), W, A); % proiectare filtru

    [H, f] = freqz(h_firpm, 1, 5000); % calc rasp frecv

    % linia 1: spectrul in coord liniare
    subplot(3, length(M_ajustat_ab), i);
    plot(f/pi, abs(H));
    hold on;
    xline(omega_p_prime/pi, '--r');
    xline(omega_s_prime/pi, '--r');
    title(sprintf('PM: M = %d', M_ajustat_ab(i)));
    xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
    grid on;

    % linia 2: spectrul in dB
    subplot(3, length(M_ajustat_ab), i + length(M_ajustat_ab));
    plot(f/pi, 20*log10(abs(H)));
    hold on;
    xline(omega_p_prime/pi, '--r');
    xline(omega_s_prime/pi, '--r');
    xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})| [dB]');
    grid on;

    % linia 3: faza
    subplot(3, length(M_ajustat_ab), i + 2*length(M_ajustat_ab));
    plot(f/pi, unwrap(angle(H)));
    xlabel('\omega / \pi'); ylabel('Fază [rad]');
    grid on;
end

% mat 3: secv pondere
matrix3 = figure('Name', 'Faza 1 b - Matricea 3: Secvențele pondere');
for i = 1:length(M_ajustat_ab)
    % proiectarea pentru {omega_p, omega_s}
    W = [0 omega_p/pi omega_s/pi 1];
    h1 = firpm(M_ajustat_ab(i), W, A);

    % proiectarea pentru {pi - omega_s, pi - omega_p}
    omega_p_prime = pi - omega_s;
    omega_s_prime = pi - omega_p;
    W_prime = [0 omega_p_prime/pi omega_s_prime/pi 1];
    h2 = firpm(M_ajustat_ab(i), W_prime, A);

    % graficul secv pondere
    subplot(length(M_ajustat_ab), 2, 2*i - 1);
    stem(h1, 'filled');
    title(sprintf('Secv. pondere PM: M = %d, {omega_p, omega_s}', M_ajustat_ab(i)));
    xlabel('n'); ylabel('h[n]');
    grid on;

    subplot(length(M_ajustat_ab), 2, 2*i);
    stem(h2, 'filled');
    title(sprintf('Secv. pondere PM: M = %d, {pi-omega_s, pi-omega_p}', M_ajustat_ab(i)));
    xlabel('n'); ylabel('h[n]');
    grid on;
end

%% c

M = PS_PRJ_2_Faza_1c(ng,ns); 

M_ajustat_cd = M - 1;

%%

% def vectorilor de frecvente și amplitudini
W = [0 omega_p/pi omega_s/pi 1];
A = [1 1 0 0];

% Proiectarea filtrelor pentru {omega_p, omega_s}
h_firls = firls(M_ajustat_cd, W, A);
h_firpm = firpm(M_ajustat_cd, W, A);

% calc răspunsului în frecv
[H_firls, f] = freqz(h_firls, 1, 5000);
[H_firpm, ~] = freqz(h_firpm, 1, 5000);

% det atenuării minime în banda de stopare
stopband = (f/pi > omega_s/pi); % Indicele benzii de stopare
atten_firls = min(20*log10(abs(H_firls(stopband))));
atten_firpm = min(20*log10(abs(H_firpm(stopband))));

% Graficul răspunsurilor în frecv
figure('Name', 'Faza 1 - Punctul c: Compararea filtrelor pentru M fixat');
plot(f/pi, 20*log10(abs(H_firls)), 'b', 'DisplayName', 'CMMP (firls)'); hold on;
plot(f/pi, 20*log10(abs(H_firpm)), 'r', 'DisplayName', 'Normă infinit (firpm)');

% Linii verticale pentru pulsatii
xline(omega_p/pi, '--k', 'DisplayName', '\omega_p');
xline(omega_s/pi, '--k', 'DisplayName', '\omega_s');

% Adăugarea titlului și legendelor
title(sprintf('Compararea filtrelor pentru M = %d: ALS = %.2f dB, APM = %.2f dB', M_ajustat_cd, atten_firls, atten_firpm));
xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})| [dB]');
grid on;
legend;

% proiectarea filtrelor pentru {pi - omega_s, pi - omega_p}
omega_p_prime = pi - omega_s;
omega_s_prime = pi - omega_p;
W_prime = [0 omega_p_prime/pi omega_s_prime/pi 1];
h_firls_prime = firls(M_ajustat_cd, W_prime, A);
h_firpm_prime = firpm(M_ajustat_cd, W_prime, A);

% calc rasp in frecvență pentru pulsatiile transformate
[H_firls_prime, f] = freqz(h_firls_prime, 1, 5000);
[H_firpm_prime, ~] = freqz(h_firpm_prime, 1, 5000);

% det atenuării minime pentru pulsatiile transformate
stopband_prime = (f/pi > omega_s_prime/pi);
atten_firls_prime = min(20*log10(abs(H_firls_prime(stopband_prime))));
atten_firpm_prime = min(20*log10(abs(H_firpm_prime(stopband_prime))));

% Graficul răspunsurilor în frecvență pentru pulsatiile transformate
figure('Name', 'Faza 1 - Punctul c: Compararea filtrelor pentru pulsatiile transformate');
plot(f/pi, 20*log10(abs(H_firls_prime)), 'b', 'DisplayName', 'CMMP (firls)'); hold on;
plot(f/pi, 20*log10(abs(H_firpm_prime)), 'r', 'DisplayName', 'Normă infinit (firpm)');

% Linii verticale pentru pulsatii transformate
xline(omega_p_prime/pi, '--k', 'DisplayName', '\pi - \omega_s');
xline(omega_s_prime/pi, '--k', 'DisplayName', '\pi - \omega_p');

% Adăugarea titlului și legendelor
title(sprintf('Compararea filtrelor pentru M = %d (pulsatii transformate): ALS = %.2f dB, APM = %.2f dB', M_ajustat_cd, atten_firls_prime, atten_firpm_prime));
xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})| [dB]');
grid on;
legend;

%% d

ws = PS_PRJ_2_Faza_1d(ng, ns);

% Definirea vectorilor pentru cele două cazuri de pulsatii
W1 = [0 omega_p/pi omega_s/pi 1];  % {ω_p, ω_s}
W2 = [0 (pi - omega_s)/pi (pi - omega_p)/pi 1];  % {π - ω_s, π - ω_p}
A = [1 1 0 0];  % Amplitudini pentru ambele cazuri

% S1: Fără ponderi
% Caz 1: {ω_p, ω_s}
h_firls_no_weight1 = firls(M_ajustat_cd, W1, A);
h_firpm_no_weight1 = firpm(M_ajustat_cd, W1, A);
[H_firls_no_weight1, f] = freqz(h_firls_no_weight1, 1, 5000);
[H_firpm_no_weight1, ~] = freqz(h_firpm_no_weight1, 1, 5000);

% Caz 2: {π - ω_s, π - ω_p}
h_firls_no_weight2 = firls(M_ajustat_cd, W2, A);
h_firpm_no_weight2 = firpm(M_ajustat_cd, W2, A);
[H_firls_no_weight2, ~] = freqz(h_firls_no_weight2, 1, 5000);
[H_firpm_no_weight2, ~] = freqz(h_firpm_no_weight2, 1, 5000);

% Calculul raporturilor dispersiilor pentru S1
dispersion_firls_no_weight1 = sum(abs(A(1) - abs(H_firls_no_weight1(f <= omega_p))).^2) + ...
                              sum(abs(A(3) - abs(H_firls_no_weight1(f >= omega_s))).^2);
dispersion_firpm_no_weight1 = sum(abs(A(1) - abs(H_firpm_no_weight1(f <= omega_p))).^2) + ...
                              sum(abs(A(3) - abs(H_firpm_no_weight1(f >= omega_s))).^2);
RLS_no_weight1 = dispersion_firls_no_weight1 / dispersion_firpm_no_weight1;

dispersion_firls_no_weight2 = sum(abs(A(1) - abs(H_firls_no_weight2(f <= pi - omega_s))).^2) + ...
                              sum(abs(A(3) - abs(H_firls_no_weight2(f >= pi - omega_p))).^2);
dispersion_firpm_no_weight2 = sum(abs(A(1) - abs(H_firpm_no_weight2(f <= pi - omega_s))).^2) + ...
                              sum(abs(A(3) - abs(H_firpm_no_weight2(f >= pi - omega_p))).^2);
RLS_no_weight2 = dispersion_firls_no_weight2 / dispersion_firpm_no_weight2;

% S2: Cu ponderi
% Caz 1: {ω_p, ω_s}
W_weights1 = [1 ws];
h_firls_weighted1 = firls(M_ajustat_cd, W1, A, W_weights1);
h_firpm_weighted1 = firpm(M_ajustat_cd, W1, A, W_weights1);
[H_firls_weighted1, ~] = freqz(h_firls_weighted1, 1, 5000);
[H_firpm_weighted1, ~] = freqz(h_firpm_weighted1, 1, 5000);

% Caz 2: {π - ω_s, π - ω_p}
W_weights2 = [1 ws];
h_firls_weighted2 = firls(M_ajustat_cd, W2, A, W_weights2);
h_firpm_weighted2 = firpm(M_ajustat_cd, W2, A, W_weights2);
[H_firls_weighted2, ~] = freqz(h_firls_weighted2, 1, 5000);
[H_firpm_weighted2, ~] = freqz(h_firpm_weighted2, 1, 5000);

% Calculul raporturilor dispersiilor pentru S2
dispersion_firls_weighted1 = sum(abs(A(1) - abs(H_firls_weighted1(f <= omega_p))).^2) + ...
                              sum(abs(A(3) - abs(H_firls_weighted1(f >= omega_s))).^2);
dispersion_firpm_weighted1 = sum(abs(A(1) - abs(H_firpm_weighted1(f <= omega_p))).^2) + ...
                              sum(abs(A(3) - abs(H_firpm_weighted1(f >= omega_s))).^2);
RLS_weighted1 = dispersion_firls_weighted1 / dispersion_firpm_weighted1;

dispersion_firls_weighted2 = sum(abs(A(1) - abs(H_firls_weighted2(f <= pi - omega_s))).^2) + ...
                              sum(abs(A(3) - abs(H_firls_weighted2(f >= pi - omega_p))).^2);
dispersion_firpm_weighted2 = sum(abs(A(1) - abs(H_firpm_weighted2(f <= pi - omega_s))).^2) + ...
                              sum(abs(A(3) - abs(H_firpm_weighted2(f >= pi - omega_p))).^2);
RLS_weighted2 = dispersion_firls_weighted2 / dispersion_firpm_weighted2;

% Grafic pentru S1: Fără ponderi
figure('Name', 'Faza 1 - Punctul d: Fără ponderi');
subplot(1, 2, 1);
plot(f/pi, abs(H_firls_no_weight1), 'b', 'DisplayName', 'CMMP (firls)'); hold on;
plot(f/pi, abs(H_firpm_no_weight1), 'r', 'DisplayName', 'Normă infinit (firpm)');
xline(omega_p/pi, '--k', 'DisplayName', '\omega_p');
xline(omega_s/pi, '--k', 'DisplayName', '\omega_s');
legend;
xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
grid on;
title(sprintf('s1: {ω_p, ω_s}, RLS = %.2f', RLS_no_weight1));

subplot(1, 2, 2);
plot(f/pi, abs(H_firls_no_weight2), 'b', 'DisplayName', 'CMMP (firls)'); hold on;
plot(f/pi, abs(H_firpm_no_weight2), 'r', 'DisplayName', 'Normă infinit (firpm)');
xline((pi - omega_s)/pi, '--k', 'DisplayName', '\pi - \omega_s');
xline((pi - omega_p)/pi, '--k', 'DisplayName', '\pi - \omega_p');
legend;
xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
grid on;
title(sprintf('s1: {π - ω_s, π - ω_p}, RLS = %.2f', RLS_no_weight2));

% Grafic pentru S2: Cu ponderi
figure('Name', 'Faza 1 - Punctul d: Cu ponderi');
subplot(1, 2, 1);
plot(f/pi, abs(H_firls_weighted1), 'b', 'DisplayName', 'CMMP (firls)'); hold on;
plot(f/pi, abs(H_firpm_weighted1), 'r', 'DisplayName', 'Normă infinit (firpm)');
xline(omega_p/pi, '--k', 'DisplayName', '\omega_p');
xline(omega_s/pi, '--k', 'DisplayName', '\omega_s');
legend;
xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
grid on;
title(sprintf('s2: {ω_p, ω_s}, RLS = %.2f', RLS_weighted1));

subplot(1, 2, 2);
plot(f/pi, abs(H_firls_weighted2), 'b', 'DisplayName', 'CMMP (firls)'); hold on;
plot(f/pi, abs(H_firpm_weighted2), 'r', 'DisplayName', 'Normă infinit (firpm)');
xline((pi - omega_s)/pi, '--k', 'DisplayName', '\pi - \omega_s');
xline((pi - omega_p)/pi, '--k', 'DisplayName', '\pi - \omega_p');
legend;
xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
grid on;
title(sprintf('s2: {π - ω_s, π - ω_p}, RLS = %.2f', RLS_weighted2));


%% e

M = 20; % Ordinul fixat al filtrului
M_ajustat_e = M - 1;

% Definirea pulsatiilor pentru FTS și FTB
omega_p_FTS = min(omega_s, omega_p); % Asigurăm ordine crescătoare
omega_s_FTS = max(omega_s, omega_p);

omega_p_FTB = omega_p;
omega_s_FTB = omega_s;
omega_p_prime_FTB = pi - omega_s;
omega_s_prime_FTB = pi - omega_p;

% Pulsatiile pentru FTB, ordonate crescător
W_FTB = sort([0 omega_p_FTB/pi omega_s_FTB/pi omega_p_prime_FTB/pi omega_s_prime_FTB/pi 1]);

% Amplitudinea pentru FTS și FTB
A_FTS = [0 0 1 1];
A_FTB = [0 0 1 1 0 0];

% Metoda 1: Proiectarea filtrelor FTS și FTB direct cu firls și firpm
% FTS
h_firls_FTS = firls(M_ajustat_e, [0 omega_p_FTS/pi omega_s_FTS/pi 1], A_FTS, 'h');
h_firpm_FTS = firpm(M_ajustat_e, [0 omega_p_FTS/pi omega_s_FTS/pi 1], A_FTS, 'h');

% FTB
h_firls_FTB = firls(M_ajustat_e, W_FTB, A_FTB, 'h');
h_firpm_FTB = firpm(M_ajustat_e, W_FTB, A_FTB, 'h');

% Metoda 2: Proiectarea filtrelor FTS și FTB pe baza lui FTJ
% FTJ
h_firls_FTJ = firls(M_ajustat_e, [0 omega_p/pi omega_s/pi 1], [1 1 0 0], 'h');
h_firpm_FTJ = firpm(M_ajustat_e, [0 omega_p/pi omega_s/pi 1], [1 1 0 0], 'h');

% FTS pe baza lui FTJ
h_firls_FTS_alt = 1 - h_firls_FTJ;
h_firpm_FTS_alt = 1 - h_firpm_FTJ;

% FTB pe baza lui FTJ
h_firls_FTB_alt = h_firls_FTJ - h_firls_FTJ;
h_firpm_FTB_alt = h_firpm_FTJ - h_firpm_FTJ;

% Calculul raspunsului în frecvență pentru toate filtrele
[H_firls_FTS, f] = freqz(h_firls_FTS, 1, 5000);
[H_firpm_FTS, ~] = freqz(h_firpm_FTS, 1, 5000);

[H_firls_FTB, ~] = freqz(h_firls_FTB, 1, 5000);
[H_firpm_FTB, ~] = freqz(h_firpm_FTB, 1, 5000);

[H_firls_FTS_alt, ~] = freqz(h_firls_FTS_alt, 1, 5000);
[H_firpm_FTS_alt, ~] = freqz(h_firpm_FTS_alt, 1, 5000);

[H_firls_FTB_alt, ~] = freqz(h_firls_FTB_alt, 1, 5000);
[H_firpm_FTB_alt, ~] = freqz(h_firpm_FTB_alt, 1, 5000);

% Norme pentru fiecare filtru
norm_firls_FTS = norm(h_firls_FTS - h_firls_FTS_alt);
norm_firpm_FTS = norm(h_firpm_FTS - h_firpm_FTS_alt);
norm_firls_FTB = norm(h_firls_FTB - h_firls_FTB_alt);
norm_firpm_FTB = norm(h_firpm_FTB - h_firpm_FTB_alt);

disp('Normele pentru secventa pondere (firls și firpm):');
disp(['Norma FTS (firls): ', num2str(norm_firls_FTS)]);
disp(['Norma FTS (firpm): ', num2str(norm_firpm_FTS)]);
disp(['Norma FTB (firls): ', num2str(norm_firls_FTB)]);
disp(['Norma FTB (firpm): ', num2str(norm_firpm_FTB)]);

% Graficul comparativ al raspunsurilor în frecvență pentru FTS
figure('Name', 'Faza 1 - Punctul e: FTS');
subplot(1, 2, 1);
plot(f/pi, abs(H_firls_FTS), 'b', 'DisplayName', 'CMMP Direct (firls)'); hold on;
plot(f/pi, abs(H_firls_FTS_alt), 'r--', 'DisplayName', 'CMMP Alternativ (1 - FTJ)');
legend;
xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
grid on;
title('FTS: CMMP Comparatie');

subplot(1, 2, 2);
plot(f/pi, abs(H_firpm_FTS), 'b', 'DisplayName', 'Norma infinita Direct (firpm)'); hold on;
plot(f/pi, abs(H_firpm_FTS_alt), 'r--', 'DisplayName', 'Norma infinita Alternativ (1 - FTJ)');
legend;
xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
grid on;
title('FTS: Norma infinita Comparatie');

% Graficul comparativ al raspunsurilor în frecvență pentru FTB
figure('Name', 'Faza 1 - Punctul e: FTB');
subplot(1, 2, 1);
plot(f/pi, abs(H_firls_FTB), 'b', 'DisplayName', 'CMMP Direct (firls)'); hold on;
plot(f/pi, abs(H_firls_FTB_alt), 'r--', 'DisplayName', 'CMMP Alternativ (FTJ - FTJ)');
legend;
xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
grid on;
title('FTB: CMMP Comparatie');

subplot(1, 2, 2);
plot(f/pi, abs(H_firpm_FTB), 'b', 'DisplayName', 'Norma infinita Direct (firpm)'); hold on;
plot(f/pi, abs(H_firpm_FTB_alt), 'r--', 'DisplayName', 'Norma infinita Alternativ (FTJ - FTJ)');
legend;
xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
grid on;
title('FTB: Norma infinita Comparatie');

% Graficul comparativ al secventelor de pondere
figure('Name', 'Faza 1 - Punctul e: Secvente de Pondere');
subplot(2, 2, 1);
stem(h_firls_FTS, 'b', 'filled'); hold on;
stem(h_firls_FTS_alt, 'r--', 'filled');
title('FTS: CMMP Comparatie');
xlabel('n'); ylabel('h[n]');
grid on;

subplot(2, 2, 2);
stem(h_firpm_FTS, 'b', 'filled'); hold on;
stem(h_firpm_FTS_alt, 'r--', 'filled');
title('FTS: Norma infinita Comparatie');
xlabel('n'); ylabel('h[n]');
grid on;

subplot(2, 2, 3);
stem(h_firls_FTB, 'b', 'filled'); hold on;
stem(h_firls_FTB_alt, 'r--', 'filled');
title('FTB: CMMP Comparatie');
xlabel('n'); ylabel('h[n]');
grid on;

subplot(2, 2, 4);
stem(h_firpm_FTB, 'b', 'filled'); hold on;
stem(h_firpm_FTB_alt, 'r--', 'filled');
title('FTB: Norma infinita Comparatie');
xlabel('n'); ylabel('h[n]');
grid on;

%% faza 2

%% a

M = 20; % Ordinul filtrului
K_values = [5, 7, 10, 13, 15]; % Întârzierea de grup

% Două perechi de pulsatii: {omega_p, omega_s} și {pi - omega_s, pi - omega_p}
omega_p_prime = pi - omega_s;
omega_s_prime = pi - omega_p;

% Matricea pentru {omega_p, omega_s}
figure('Name', 'Faza 2 a - Matricea 1: {omega_p, omega_s}');
for i = 1:length(K_values)
    K = K_values(i);

    % Proiectarea filtrului
    h = firls_FTJ_c_modif(M, omega_p/pi, omega_s/pi, K);

    % Răspuns în frecvență
    [H, f] = freqz(h, 1, 5000, 'whole');

    % Întârzierea de grup
    [grp_delay, ~] = grpdelay(h, 1, 5000, 'whole');

    % Graficul caracteristicii de frecvență
    subplot(3, length(K_values), i);
    plot(f - pi, abs(fftshift(H)));
    hold on;
    xline(omega_p - pi, '--r');
    xline(omega_s - pi, '--r');
    title(sprintf('K = %d', K));
    xlabel('\omega'); ylabel('|H(e^{j\omega})|');
    grid on;

    % Graficul întârzierii de grup
    subplot(3, length(K_values), i + length(K_values));
    plot(f - pi, fftshift(grp_delay));
    hold on;
    xline(omega_p - pi, '--r');
    xline(omega_s - pi, '--r');
    xlabel('\omega'); ylabel('Întârziere de grup');
    grid on;

    % Graficul secvenței pondere
    subplot(3, length(K_values), i + 2 * length(K_values));
    stem(h, 'filled');
    xlabel('n'); ylabel('h[n]');
    grid on;
end

% Matricea pentru {pi - omega_s, pi - omega_p}
figure('Name', 'Faza 2 a - Matricea 2: {pi - omega_s, pi - omega_p}');
for i = 1:length(K_values)
    K = K_values(i);

    % Proiectarea filtrului
    h = firls_FTJ_c_modif(M, omega_p_prime/pi, omega_s_prime/pi, K);

    % Răspuns în frecvență
    [H, f] = freqz(h, 1, 5000, 'whole');

    % Întârzierea de grup
    [grp_delay, ~] = grpdelay(h, 1, 5000, 'whole');

    % Graficul caracteristicii de frecvență
    subplot(3, length(K_values), i);
    plot(f - pi, abs(fftshift(H)));
    hold on;
    xline(omega_p_prime - pi, '--r');
    xline(omega_s_prime - pi, '--r');
    title(sprintf('K = %d', K));
    xlabel('\omega'); ylabel('|H(e^{j\omega})|');
    grid on;

    % Graficul întârzierii de grup
    subplot(3, length(K_values), i + length(K_values));
    plot(f - pi, fftshift(grp_delay));
    hold on;
    xline(omega_p_prime - pi, '--r');
    xline(omega_s_prime - pi, '--r');
    xlabel('\omega'); ylabel('Întârziere de grup');
    grid on;

    % Graficul secvenței pondere
    subplot(3, length(K_values), i + 2 * length(K_values));
    stem(h, 'filled');
    xlabel('n'); ylabel('h[n]');
    grid on;
end

%% b

% Intervalul pentru ordinul filtrului
M_values = M_ajustat_ab;
K_values = M_values / 2; % Întârzierea de grup K este M/2

% Matricea 1: Caracteristici spectrale pentru {omega_p, omega_s}
figure('Name', 'Faza 2 b - Matricea 1: Caracteristici spectrale');
for i = 1:length(M_values)
    M = M_values(i);
    K = K_values(i);

    % Proiectarea filtrului
    [h, pr] = firls_FTJ_c_modif(M, omega_p/pi, omega_s/pi, K);

    % Calculul răspunsului în frecvență
    [H, f] = freqz(h, 1, 5000, 'whole');

    % Graficul caracteristicii de frecvență
    subplot(3, length(M_values), i);
    plot(f - pi, abs(fftshift(H)));
    hold on;
    xline(omega_p - pi, '--r');
    xline(omega_s - pi, '--r');
    title(sprintf('M = %d, PR = %.2f%%', M, pr));
    xlabel('\omega'); ylabel('|H(e^{j\omega})|');
    grid on;

    % Graficul întârzierii de grup
    [grp_delay, ~] = grpdelay(h, 1, 5000, 'whole');
    subplot(3, length(M_values), i + length(M_values));
    plot(f - pi, fftshift(grp_delay));
    hold on;
    xline(omega_p - pi, '--r');
    xline(omega_s - pi, '--r');
    xlabel('\omega'); ylabel('Întârziere de grup');
    grid on;

    % Graficul secvenței pondere
    subplot(3, length(M_values), i + 2 * length(M_values));
    stem(h, 'filled');
    xlabel('n'); ylabel('h[n]');
    grid on;
end

% Matricea 2: Caracteristici spectrale pentru {pi - omega_s, pi - omega_p}
omega_p_prime = pi - omega_s;
omega_s_prime = pi - omega_p;

figure('Name', 'Faza 2 b - Matricea 2: Caracteristici spectrale {pi - omega_s, pi - omega_p}');
for i = 1:length(M_values)
    M = M_values(i);
    K = K_values(i);

    % Proiectarea filtrului
    [h, pr] = firls_FTJ_c_modif(M, omega_p_prime/pi, omega_s_prime/pi, K);

    % Calculul răspunsului în frecvență
    [H, f] = freqz(h, 1, 5000, 'whole');

    % Graficul caracteristicii de frecvență
    subplot(3, length(M_values), i);
    plot(f - pi, abs(fftshift(H)));
    hold on;
    xline(omega_p_prime - pi, '--r');
    xline(omega_s_prime - pi, '--r');
    title(sprintf('M = %d, PR = %.2f%%', M, pr));
    xlabel('\omega'); ylabel('|H(e^{j\omega})|');
    grid on;

    % Graficul întârzierii de grup
    [grp_delay, ~] = grpdelay(h, 1, 5000, 'whole');
    subplot(3, length(M_values), i + length(M_values));
    plot(f - pi, fftshift(grp_delay));
    hold on;
    xline(omega_p_prime - pi, '--r');
    xline(omega_s_prime - pi, '--r');
    xlabel('\omega'); ylabel('Întârziere de grup');
    grid on;

    % Graficul secvenței pondere
    subplot(3, length(M_values), i + 2 * length(M_values));
    stem(h, 'filled');
    xlabel('n'); ylabel('h[n]');
    grid on;
end

%% c

M = 20; % Ordinul fixat al filtrului
M_ajustat = M - 1;
K = M / 2; % Întârzierea de grup pentru FTS și FTB

% Pulsatiile pentru FTS
omega_p_FTS = min(omega_s, omega_p); % Banda de stopare devine banda de trecere
omega_s_FTS = max(omega_s, omega_p);

% Pulsatiile pentru FTB
omega_p_FTB = omega_p;
omega_s_FTB = omega_s;
omega_p_prime_FTB = pi - omega_s;
omega_s_prime_FTB = pi - omega_p;

% Vectorul de frecvențe pentru FTB
W_FTB = sort([0 omega_p_FTB/pi omega_s_FTB/pi omega_p_prime_FTB/pi omega_s_prime_FTB/pi 1]);

% Amplitudinile pentru FTS și FTB
A_FTS = [0 0 1 1];
A_FTB = [0 0 1 1 0 0];

% Proiectarea filtrelor
% FTS
[h_FTS, pr_FTS] = firls_FTS_c(M_ajustat, omega_p_FTS/pi, omega_s_FTS/pi, K);

% FTB
[h_FTB, pr_FTB] = firls_FTB_c(M_ajustat, W_FTB, K);

% Calculul răspunsurilor în frecvență
[H_FTS, f] = freqz(h_FTS, 1, 5000, 'whole');
[H_FTB, ~] = freqz(h_FTB, 1, 5000, 'whole');

% Graficul răspunsurilor în frecvență
figure('Name', 'Faza 2 - Punctul c: Răspunsuri în frecvență');
subplot(2, 1, 1);
plot(f/pi, abs(H_FTS), 'b', 'DisplayName', 'FTS'); hold on;
xline(omega_p_FTS/pi, '--r', 'DisplayName', '\omega_p');
xline(omega_s_FTS/pi, '--g', 'DisplayName', '\omega_s');
legend;
xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
title('Răspunsul în frecvență: FTS');
grid on;

subplot(2, 1, 2);
plot(f/pi, abs(H_FTB), 'r', 'DisplayName', 'FTB'); hold on;
xline(omega_p_FTB/pi, '--r', 'DisplayName', '\omega_p');
xline(omega_s_FTB/pi, '--g', 'DisplayName', '\omega_s');
xline(omega_p_prime_FTB/pi, '--c', 'DisplayName', '\pi - \omega_s');
xline(omega_s_prime_FTB/pi, '--m', 'DisplayName', '\pi - \omega_p');
legend;
xlabel('\omega / \pi'); ylabel('|H(e^{j\omega})|');
title('Răspunsul în frecvență: FTB');
grid on;

% Graficul secvențelor de pondere
figure('Name', 'Faza 2 - Punctul c: Secvențe de pondere');
subplot(2, 1, 1);
stem(h_FTS, 'filled');
title('Secvența de pondere: FTS');
xlabel('n'); ylabel('h[n]');
grid on;

subplot(2, 1, 2);
stem(h_FTB, 'filled');
title('Secvența de pondere: FTB');
xlabel('n'); ylabel('h[n]');
grid on;
% Afișarea performanțelor relative
fprintf('Performanța relativă pentru FTS: %.2f%%\n', pr_FTS);
fprintf('Performanța relativă pentru FTB: %.2f%%\n', pr_FTB);