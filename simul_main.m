%% Unified simulation: Binomial LRT vs Two-way ANOVA vs Normal-ILRT
clear; clc;
rng(123);

%% ---------------- Basic settings ----------------
alpha      = 0.05;          % significance level
num_simul  = 1000;          % number of simulation runs
num_seeds  = 30;            % n_ijk (number of seeds per replicate)
rep_num    = 4;             % N_ij (number of replicates per cell)

%% ---------------- Scenario (p_ij) ----------------
% Choose one of the following rate_vec scenarios.

%  scenario 1 - many zero cells
%rate_vec = [0,0,0.2,0.5,  0,0,0.2,0.5,  0,0,0.2,0.5,  0,0,0.2,0.5];
rate_vec = [0,0,0.2,0.3,  0,0,0.2,0.3,  0,0,0.2,0.3,  0,0,0.2,0.3];

% % scenario 1 - no zero cells
%rate_vec = [0.1,0.2,0.4,0.6,  0.1,0.2,0.4,0.6,  0.1,0.2,0.4,0.6,  0.1,0.2,0.4,0.6];

% % scenario 2-1
% rate_vec = [0,0.1,0.3,0.3,  0.2,0.3,0.5,0.5,  0.5,0.6,0.8,0.8,  0.7,0.8,1,1];

% % scenario 2-2
% rate_vec = [0.2,0.3,0.5,0.5,  0.2,0.3,0.5,0.5,  0.4,0.5,0.7,0.7, 0.5,0.6,0.8,0.8];

% scenario 3 (default)
%%rate_vec = [0,0.1,0.4,0.4,  0.2,0.3,0.6,0.6,  0.5,0.6,0.6,0.7,  0.6,0.7,1,1];

P = reshape(rate_vec, 4, 4);        % rows = Factor A, cols = Factor B
[I, J]  = size(P);
IJ       = I*J;
n_vec    = rep_num * ones(1, IJ);   % number of replicates per cell
N_total  = sum(n_vec);              % total number of replicates

%% ---------------- Storage for results ----------------
% 1) Binomial LRT (ApB vs AB, A vs ApB, B vs ApB)
log_diff_ApBvsAB = zeros(num_simul,1);  % Interaction
log_diff_AvsApB  = zeros(num_simul,1);  % B effect (ApB vs A)
log_diff_BvsApB  = zeros(num_simul,1);  % A effect (ApB vs B)

% 2) Two-way ANOVA
Pval_A  = zeros(num_simul,1);
Pval_B  = zeros(num_simul,1);
Pval_AB = zeros(num_simul,1);

% 3) Normal-ILRT
T_A     = nan(num_simul,1);
T_B     = nan(num_simul,1);
T_Int   = nan(num_simul,1);

% Critical values (chi-square)
crit_A   = chi2inv(1 - alpha, I-1);
crit_B   = chi2inv(1 - alpha, J-1);
crit_Int = chi2inv(1 - alpha, (I-1)*(J-1));

% Small-sample adjustment (adf): per the original paper
if ((I + J) > 9)
    adf = ones(1, IJ);
    for nn = 1:IJ
        if n_vec(nn) == 5, adf(nn) = 0.947 + 0.00234*n_vec(nn); end
    end
else
    adf = ones(1, IJ);
end

%% ---------------- Simulation loop ----------------
for iter = 1:num_simul
    fprintf('Simulation %d / %d\n', iter, num_simul);

    %% Data generation (table form)
    % raw_data: Fac_A, Fac_B, g_ijk, n_ijk (one row per replicate)
    raw_data = table([], [], [], [], 'VariableNames', {'Fac_A','Fac_B','g_ijk','n_ijk'});
    for i = 1:I
        for j = 1:J
            p_ij = P(i,j);
            for r = 1:rep_num
                g = sum(binornd(1, p_ij, [1, num_seeds])); % Binomial(num_seeds, p_ij)
                new_row = table(i, j, g, num_seeds, ...
                    'VariableNames', {'Fac_A','Fac_B','g_ijk','n_ijk'});
                raw_data = [raw_data; new_row]; %#ok<AGROW>
            end
        end
    end
    df = raw_data;

    %% ===== 1) Binomial LRT =====
    % M_{AB}: cellwise MLEs
    p_hat_AB = varfun(@(x) sum(x), df, 'InputVariables', {'g_ijk','n_ijk'}, ...
        'GroupingVariables', {'Fac_A','Fac_B'});
    p_hat_AB.p_hat_ij = p_hat_AB.Fun_g_ijk ./ p_hat_AB.Fun_n_ijk;

    df_AB = innerjoin(df, p_hat_AB(:, {'Fac_A','Fac_B','p_hat_ij'}), 'Keys', {'Fac_A','Fac_B'});
    df_AB.loglik = log(binopdf(df_AB.g_ijk, df_AB.n_ijk, df_AB.p_hat_ij));
    total_LL_AB  = sum(df_AB.loglik);

    % M_A: row-wise (A-only) MLEs
    p_hat_A = varfun(@(x) sum(x), df, 'InputVariables', {'g_ijk','n_ijk'}, ...
        'GroupingVariables', {'Fac_A'});
    p_hat_A.p_hat_ij = p_hat_A.Fun_g_ijk ./ p_hat_A.Fun_n_ijk;

    df_A = innerjoin(df, p_hat_A(:, {'Fac_A','p_hat_ij'}), 'Keys', {'Fac_A'});
    df_A.loglik = log(binopdf(df_A.g_ijk, df_A.n_ijk, df_A.p_hat_ij));
    total_LL_A  = sum(df_A.loglik);

    % M_B: column-wise (B-only) MLEs
    p_hat_B = varfun(@(x) sum(x), df, 'InputVariables', {'g_ijk','n_ijk'}, ...
        'GroupingVariables', {'Fac_B'});
    p_hat_B.p_hat_ij = p_hat_B.Fun_g_ijk ./ p_hat_B.Fun_n_ijk;

    df_B = innerjoin(df, p_hat_B(:, {'Fac_B','p_hat_ij'}), 'Keys', {'Fac_B'});
    df_B.loglik = log(binopdf(df_B.g_ijk, df_B.n_ijk, df_B.p_hat_ij));
    total_LL_B  = sum(df_B.loglik);

    % M_{A+B}: additive constraints (fmincon)
    % (Keep original structure; avoid variable name collisions)
    Aeq4 = eye(I); Aeq4(2:I) = 1;
    L1   = kron(Aeq4,             ones(J,1));
    L2   = kron(ones(I,1), [zeros(1,J-1); eye(J-1)]);
    L    = [L1, L2];
    constr_L = L(2:end, :);

    Acon = [constr_L; -constr_L];
    Bcon = [ones(size(constr_L,1),1); zeros(size(constr_L,1),1)];

    % Initial value: based on model B
    tmpB       = p_hat_B.p_hat_ij;
    theta_ini  = [tmpB(1); zeros(I-1,1); tmpB(2:J) - tmpB(1)];

    lb = [0;  -Inf(I + J - 2, 1)];
    ub = [1;   Inf(I + J - 2, 1)];

    eval_f0 = @(theta) nlogApB_fun(theta, df, L, rep_num);
    opts_fmincon = optimoptions('fmincon','Algorithm','interior-point','Display','off');
    theta = fmincon(eval_f0, theta_ini, Acon, Bcon, [], [], lb, ub, [], opts_fmincon);

    p_hat_ApB       = kron(L * theta, ones(rep_num,1));
    df_ApB          = df;
    df_ApB.p_hat_ij = p_hat_ApB;
    df_ApB.loglik   = log(binopdf(df_ApB.g_ijk, df_ApB.n_ijk, p_hat_ApB));
    total_LL_ApB    = sum(df_ApB.loglik);

    % LRT statistics (correct mapping)
    %  - A effect:      M_{A+B} vs M_B      ¡æ log_diff_BvsApB, df = I-1
    %  - B effect:      M_{A+B} vs M_A      ¡æ log_diff_AvsApB, df = J-1
    %  - Interaction:   M_{AB}  vs M_{A+B}  ¡æ log_diff_ApBvsAB, df = (I-1)(J-1)
    log_diff_ApBvsAB(iter) = 2 * (total_LL_AB  - total_LL_ApB);
    log_diff_AvsApB(iter)  = 2 * (total_LL_ApB - total_LL_A);
    log_diff_BvsApB(iter)  = 2 * (total_LL_ApB - total_LL_B);

    %% ===== 2) Two-way ANOVA =====
    data_ANOVA = reshape(df.g_ijk ./ df.n_ijk, [J*rep_num, I]);
    [pval_model, ~, ~] = anova2(data_ANOVA, rep_num, 'off');
    Pval_A(iter)  = pval_model(1);
    Pval_B(iter)  = pval_model(2);
    Pval_AB(iter) = pval_model(3);

    %% ===== 3) Normal-ILRT (f3/f4/fA0/fB0) =====
    % Amat: each column = cell (t), each row = replicate (Y = g/n)
    Amat = NaN(rep_num, IJ);
    mu_hat_cell = zeros(1, IJ);
    t = 1;
    for i = 1:I
        for j = 1:J
            idx = find(df.Fac_A==i & df.Fac_B==j);
            y = df.g_ijk(idx) ./ df.n_ijk(idx);       % length = rep_num
            Amat(1:rep_num, t) = y(:);
            mu_hat_cell(t)     = mean(y);
            t = t + 1;
        end
    end

    mu_hat0 = sum(sum(Amat, 'omitnan'), 'omitnan') / N_total;

    AA = reshape(sum(Amat, 'omitnan'), J, I).';  % I¡¿J: cell sums
    BB = reshape(n_vec,                 J, I).'; % I¡¿J: cell sizes

    row_sum = sum(AA, 2);  row_n = sum(BB, 2);
    col_sum = sum(AA, 1);  col_n = sum(BB, 1);

    a0 = (row_sum ./ row_n).';      % 1¡¿I
    b0 =  col_sum ./ col_n;         % 1¡¿J

    opts_fmin = optimset('Display','off');

    % f3: no interaction (A and B only)
    ab0_f3 = [mu_hat0, a0 - mu_hat0, b0 - mu_hat0];
    [b3, F3VAL, FLAG3] = fminsearch(@(z)(-f3(n_vec, Amat, I, J, adf, z)), ab0_f3, opts_fmin);
    if FLAG3 < 1, continue; end

    % f4: with interaction
    gm0 = zeros(1, IJ);
    ab1_f4 = [b3, gm0];
    [~, F4VAL, FLAG4] = fminsearch(@(z)(-f4(n_vec, Amat, I, J, adf, z)), ab1_f4, opts_fmin);

    if FLAG4 < 1
        % If it fails, retry with residual-based initialization
        aa = b3(2:I+1);
        bb = b3(I+2:I+J+1);
        t = 1; gm_alt = zeros(1, IJ);
        for i = 1:I
            for j = 1:J
                gm_alt(t) = mu_hat_cell(t) - b3(1) - aa(i) - bb(j);
                t = t + 1;
            end
        end
        ab1_f4 = [b3, gm_alt];
        [~, F4VAL, FLAG4] = fminsearch(@(z)(-f4(n_vec, Amat, I, J, adf, z)), ab1_f4, opts_fmin);
        if FLAG4 < 1, continue; end
    end

    % fA0 (no A): [mu, beta(1:J)] vs f3
    ab0_A = [mu_hat0, b0 - mu_hat0];
    [~, F0A, FLAG0A] = fminsearch(@(z)(-fA0(n_vec, Amat, I, J, adf, z)), ab0_A, opts_fmin);
    if FLAG0A < 1, continue; end

    % fB0 (no B): [mu, alpha(1:I)] vs f3
    ab0_B = [mu_hat0, a0 - mu_hat0];
    [~, F0B, FLAG0B] = fminsearch(@(z)(-fB0(n_vec, Amat, I, J, adf, z)), ab0_B, opts_fmin);
    if FLAG0B < 1, continue; end

    % ILRT statistics
    T_Int(iter) = -2 * (-F3VAL + F4VAL);   % Interaction: H0=f3, H1=f4
    T_A(iter)   = -2 * (-F0A    + F3VAL);  % A:          H0=fA0, H1=f3
    T_B(iter)   = -2 * (-F0B    + F3VAL);  % B:          H0=fB0, H1=f3
end

%% ---------------- Detection rate summary ----------------
% 1) Binomial LRT 
det_A_binLRT   = mean(log_diff_BvsApB  > crit_A);     % A effect, df=I-1
det_B_binLRT   = mean(log_diff_AvsApB  > crit_B);     % B effect, df=J-1
det_Int_binLRT = mean(log_diff_ApBvsAB > crit_Int);   % Interaction, df=(I-1)(J-1)

% 2) Two-way ANOVA
det_A_anova    = mean(Pval_A  < alpha);
det_B_anova    = mean(Pval_B  < alpha);
det_Int_anova  = mean(Pval_AB < alpha);

% 3) Normal-ILRT
det_A_ILRT     = mean(T_A   > crit_A,   'omitnan');
det_B_ILRT     = mean(T_B   > crit_B,   'omitnan');
det_Int_ILRT   = mean(T_Int > crit_Int, 'omitnan');

fprintf('\n=== Detection rates (alpha=%.2f, simul=%d) ===\n', alpha, num_simul);
fprintf('Binomial LRT  :  A=%.4f  B=%.4f  Int=%.4f\n', det_A_binLRT, det_B_binLRT, det_Int_binLRT);
fprintf('Two-way ANOVA :  A=%.4f  B=%.4f  Int=%.4f\n', det_A_anova,  det_B_anova,  det_Int_anova);
fprintf('Normal-ILRT   :  A=%.4f  B=%.4f  Int=%.4f\n', det_A_ILRT,   det_B_ILRT,   det_Int_ILRT);
