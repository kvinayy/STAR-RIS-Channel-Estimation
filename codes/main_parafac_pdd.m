%%=======================================================
%  Combined PARAFAC Channel Estimation + PDD Optimization
%  Integrates Paper 1 (PARAFAC) with Paper 2 (PDD-BCD)
%%========================================================
clc
clear all
close all

% Add paths (adjust to your folder structure)
addpath("functions\");  % For PDD functions

%% ============ SETUP CVX SOLVER ============
% Choose one:
% cvx_solver mosek     % Use this if MOSEK is installed
cvx_solver sdpt3   % Or use this free alternative

%% ============ PARAFAC PARAMETERS (Paper 1) ============
K1 = 4;
K2 = 4;
N1 = 4;
N2 = 4;
K = K1 * K2;        % BS antennas (16)
N = N1 * N2;        % RIS elements (16)

M = 6;              % Reflection-side users
M1 = 6;             % Transmission-side users

T = 30;             % Time slots
P = 16;             % Pilot patterns
var_channel = 1;    % Channel variance
iter = 30;          % PARAFAC iterations

% Choose SNR for channel estimation
SNR_dB = 20;        % You can change this value
var_noise = 10^(-0.1*SNR_dB);

fprintf('========================================\n');
fprintf('STEP 1: PARAFAC Channel Estimation\n');
fprintf('========================================\n');
fprintf('SNR = %d dB\n', SNR_dB);

%% ============ GENERATE TRUE CHANNELS ============
G_true  = sqrt(var_channel/2) * (randn(K,N) + 1i*randn(K,N));
Hr_true = sqrt(var_channel/2) * (randn(N,M) + 1i*randn(N,M));
Ht_true = sqrt(var_channel/2) * (randn(N,M1) + 1i*randn(N,M1));

% Normalization (remove scaling ambiguity)
Hr_true(:,1) = 1;
Ht_true(:,1) = 1;

%% ============ GENERATE PILOT SIGNALS ============
[X, X_inv] = Transceiver(M, T);
[X1, X1_inv] = Transceiver1(M1, T);

%% ============ GENERATE PILOT PHASES ============
[Phi] = Phase_Generate(P, N);
[Theta] = Phase_Generate1(P, N);

%% ============ SIMULATE RECEIVED SIGNALS ============
% Reflection side
rec_y_TEMP = zeros(K, M, P);
for p = 1:P
    noise = sqrt(var_noise/2) * (randn(K,T) + 1i*randn(K,T));
    rec_y = G_true * diag(Phi(p,:)) * Hr_true * X + noise;
    rec_y_TEMP(:,:,p) = rec_y * X_inv;
end

% Transmission side
rec_y_t_TEMP = zeros(K, M1, P);
for p = 1:P
    noise = sqrt(var_noise/2) * (randn(K,T) + 1i*randn(K,T));
    rec_y_t = G_true * diag(Theta(p,:)) * Ht_true * X1 + noise;
    rec_y_t_TEMP(:,:,p) = rec_y_t * X1_inv;
end

%% ============ TENSOR UNFOLDING ============
% For Hr estimation (Mode-2 unfolding)
Z_KP_M = zeros(K*P, M);
for m = 1:M
    for p = 1:P
        for k = 1:K
            Z_KP_M((p-1)*K+k, m) = rec_y_TEMP(k,m,p);
        end
    end
end

% For Ht estimation
Z_KP_M1 = zeros(K*P, M1);
for m1 = 1:M1
    for p1 = 1:P
        for k1 = 1:K
            Z_KP_M1((p1-1)*K+k1, m1) = rec_y_t_TEMP(k1,m1,p1);
        end
    end
end

% For G estimation (Mode-1 unfolding)
Z_PM_K = zeros(P*M, K);
for m = 1:M
    for p = 1:P
        for k = 1:K
            Z_PM_K((m-1)*P+p, k) = rec_y_TEMP(k,m,p);
        end
    end
end

%% ============ PARAFAC ITERATIVE ESTIMATION ============
Hr_est = zeros(N, M, iter+1);
Ht_est = zeros(N, M1, iter+1);
G_est = zeros(K, N, iter+1);

% Random initialization
G_est(:,:,1) = sqrt(var_channel/2) * (randn(K,N) + 1i*randn(K,N));

fprintf('Running PARAFAC iterations...\n');
for i = 2:iter
    % Estimate Hr
    A2 = kr(Phi, G_est(:,:,i-1));
    A2_inv = inv(A2'*A2) * A2';
    Hr_est(:,:,i) = A2_inv * Z_KP_M;
    
    % Remove scaling ambiguity
    for n = 1:N
        Hr_est(n,:,i) = Hr_est(n,:,i) / Hr_est(n,1,i);
    end
    
    % Estimate Ht
    Ht_est(:,:,i) = A2_inv * Z_KP_M1;
    
    % Remove scaling ambiguity
    for n = 1:N
        Ht_est(n,:,i) = Ht_est(n,:,i) / Ht_est(n,1,i);
    end
    
    % Estimate G
    A1 = kr(Hr_est(:,:,i).', Phi);
    A1_inv = inv(A1'*A1) * A1';
    G_est(:,:,i) = (A1_inv * Z_PM_K).';
    
    % Check convergence
    fit = norm(Z_KP_M - kr(Phi, G_est(:,:,i)) * Hr_est(:,:,i), 'fro')^2;
    if i > 2
        delta = abs((fit_prev - fit) / fit);
        if delta < 1e-5
            break;
        end
    end
    fit_prev = fit;
end

% Extract final estimates
G_final = G_est(:,:,i);
Hr_final = Hr_est(:,:,i);
Ht_final = Ht_est(:,:,i);

% Compute NMSE
nmse_G = norm(G_true - G_final, 'fro')^2 / norm(G_true, 'fro')^2;
nmse_Hr = norm(Hr_true - Hr_final, 'fro')^2 / norm(Hr_true, 'fro')^2;
nmse_Ht = norm(Ht_true - Ht_final, 'fro')^2 / norm(Ht_true, 'fro')^2;

fprintf('PARAFAC converged at iteration %d\n', i);
fprintf('NMSE - G: %.4f, Hr: %.4f, Ht: %.4f\n', nmse_G, nmse_Hr, nmse_Ht);

%% ============ PREPARE CHANNELS FOR PDD OPTIMIZER ============
fprintf('\n========================================\n');
fprintf('STEP 2: Preparing Channels for PDD\n');
fprintf('========================================\n');

% Initialize PDD parameters
para_pdd = para_init();

% *** CRITICAL: Adjust parameters to match PARAFAC dimensions ***
para_pdd.M = K;                    % BS antennas (16 from PARAFAC)
para_pdd.N = N;                    % RIS elements (16 from PARAFAC)
para_pdd.K = M + M1;               % Total users (12 = 6 reflection + 6 transmission)

% *** CRITICAL: Channel format conversion ***
% PARAFAC: G is K×N (BS × RIS)
% PDD expects: G is N×M (RIS × BS)
G_pdd = G_final.';  % TRANSPOSE!

% Combine Hr and Ht into single h matrix
% PDD expects: h is N×K (RIS × all users)
h_pdd = [Hr_final, Ht_final];  % Concatenate: N×(M+M1)

fprintf('Channel dimensions for PDD:\n');
fprintf('  G_pdd: %d × %d (RIS × BS)\n', size(G_pdd,1), size(G_pdd,2));
fprintf('  h_pdd: %d × %d (RIS × Users)\n', size(h_pdd,1), size(h_pdd,2));

%% ============ RUN PDD OPTIMIZATION ============
fprintf('\n========================================\n');
fprintf('STEP 3: Running PDD Optimization\n');
fprintf('========================================\n');

[rate, W_opt, theta_t_opt, theta_r_opt, phas_diff_all, sum_rate_all] = ...
    algorithm_PDD(para_pdd, G_pdd, h_pdd);

%% ============ EXTRACT OPTIMAL PARAMETERS ============
fprintf('\n========================================\n');
fprintf('RESULTS: Optimal STAR-RIS Configuration\n');
fprintf('========================================\n');

% Extract amplitudes from phase-shift coefficients
beta_t_opt = abs(theta_t_opt);
beta_r_opt = abs(theta_r_opt);

% Extract phases
phase_t_opt = angle(theta_t_opt);
phase_r_opt = angle(theta_r_opt);

fprintf('Optimal Beamforming Matrix W: %d × %d\n', size(W_opt,1), size(W_opt,2));
fprintf('Optimal Transmission Phases θ_t: %d elements\n', length(phase_t_opt));
fprintf('Optimal Reflection Phases θ_r: %d elements\n', length(phase_r_opt));
fprintf('Optimal Transmission Amplitudes β_t: %d elements\n', length(beta_t_opt));
fprintf('Optimal Reflection Amplitudes β_r: %d elements\n', length(beta_r_opt));

fprintf('\nUser Rates (bits/s/Hz):\n');
for k = 1:para_pdd.K
    fprintf('  User %d: %.4f\n', k, rate(k));
end
fprintf('Total Sum Rate: %.4f bits/s/Hz\n', sum(rate));

%% ============ VISUALIZE RESULTS ============
figure('Position', [100 100 1200 800]);

% Subplot 1: PDD Convergence
subplot(2,2,1);
plot(1:length(sum_rate_all), sum_rate_all, '-bo', 'LineWidth', 1.5);
xlabel('Cumulative BCD Iterations');
ylabel('Sum Rate (bit/s/Hz)');
title('PDD Algorithm Convergence');
grid on;

% Subplot 2: Phase Difference Evolution
subplot(2,2,2);
plot(phas_diff_all', 'LineWidth', 1);
xlabel('Cumulative BCD Iterations');
ylabel('Phase Shift Difference (rad)');
title('Phase Difference θ_r - θ_t');
ylim([0, 2*pi]);
yticks([0, pi/2, pi, 3*pi/2, 2*pi]);
yticklabels({'0', 'π/2', 'π', '3π/2', '2π'});
grid on;

% Subplot 3: Optimal Amplitudes
subplot(2,2,3);
bar(1:N, [beta_t_opt, beta_r_opt]);
xlabel('RIS Element Index');
ylabel('Amplitude');
legend({'β_t (Transmission)', 'β_r (Reflection)'});
title('Optimal STAR-RIS Amplitudes');
grid on;

% Subplot 4: Optimal Phases
subplot(2,2,4);
plot(1:N, phase_t_opt, '-o', 'LineWidth', 1.5); hold on;
plot(1:N, phase_r_opt, '-s', 'LineWidth', 1.5);
xlabel('RIS Element Index');
ylabel('Phase (rad)');
legend({'θ_t (Transmission)', 'θ_r (Reflection)'});
title('Optimal STAR-RIS Phase Shifts');
grid on;
ylim([-pi, pi]);
yticks([-pi, -pi/2, 0, pi/2, pi]);
yticklabels({'-π', '-π/2', '0', 'π/2', 'π'});

sgtitle(['Combined PARAFAC + PDD Results (SNR = ' num2str(SNR_dB) ' dB)']);

%% ============ SAVE RESULTS ============
results.G_estimated = G_final;
results.Hr_estimated = Hr_final;
results.Ht_estimated = Ht_final;
results.nmse_G = nmse_G;
results.nmse_Hr = nmse_Hr;
results.nmse_Ht = nmse_Ht;
results.W_optimal = W_opt;
results.theta_t_optimal = theta_t_opt;
results.theta_r_optimal = theta_r_opt;
results.beta_t_optimal = beta_t_opt;
results.beta_r_optimal = beta_r_opt;
results.user_rates = rate;
results.sum_rate = sum(rate);
results.SNR_dB = SNR_dB;

save('combined_results.mat', 'results');
fprintf('\nResults saved to combined_results.mat\n');