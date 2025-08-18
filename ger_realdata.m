%% Real-data analysis for two-factor germination: Binomial LRT
% This script computes LRT p-values for:
%   - Interaction:   M_{AB}  vs M_{A+B}
%   - Factor A main: M_{A+B} vs M_{B}
%   - Factor B main: M_{A+B} vs M_{A}
%
% Data layout (df): one row per replicate with variables
%   Fac_A, Fac_B  : factor levels (positive integers)
%   g_ijk         : number of germinated seeds in the replicate
%   n_ijk         : number of tested seeds in the replicate (binomial trials)
%
% Requirements:
%   - The function nlogApB_fun_un.m implements the negative log-likelihood
%     for the additive (A+B) model with possibly unequal replicates per cell.

clear;
df = readtable('sheepgrass.csv');

Fac_A = table2array(df(:, 'Fac_A'));
Fac_B = table2array(df(:, 'Fac_B'));

%% Count the number of replicates per (A,B) cell (unbalanced allowed)
% Unique (Fac_A, Fac_B) combinations present in the data
unique_combinations = unique([Fac_A, Fac_B], 'rows');

% Number of rows (replicates) per unique (A,B) combination
counts = zeros(size(unique_combinations, 1), 1);
for i = 1:size(unique_combinations, 1)
    counts(i) = sum(Fac_A == unique_combinations(i,1) & Fac_B == unique_combinations(i,2));
end

% Summary table: replicate counts per cell
ctable = table(unique_combinations(:,1), unique_combinations(:,2), counts, ...
                     'VariableNames', {'Fac_A', 'Fac_B', 'Count'});

rep_vec = ctable.Count;            % vector of replicate counts per (A,B) cell
I = length(unique(Fac_A));         % number of levels for Factor A
J = length(unique(Fac_B));         % number of levels for Factor B

%% Global null model M_0: common probability p across all observations
% MLE under M_0 is the overall proportion of successes
p0 = sum(df.g_ijk) / sum(df.n_ijk);
p_hat_ij = p0 * ones(sum(rep_vec), 1);

df_0 = df;
df_0.p_hat_ij = p_hat_ij;
df_0.log_likelihood_ijk = log(binopdf(df_0.g_ijk, df_0.n_ijk, p_hat_ij));
total_log_likelihood_0 = sum(df_0.log_likelihood_ijk); 

%% Saturated cell model M_{AB}: one probability per (A,B) cell
p_hat_AB = varfun(@(x) sum(x), df, 'InputVariables', {'g_ijk', 'n_ijk'}, ...
    'GroupingVariables', {'Fac_A', 'Fac_B'});
p_hat_AB.p_hat_ij = p_hat_AB.Fun_g_ijk ./ p_hat_AB.Fun_n_ijk;

df_AB = innerjoin(df, p_hat_AB(:, {'Fac_A', 'Fac_B', 'p_hat_ij'}), 'Keys', {'Fac_A', 'Fac_B'});
df_AB.log_likelihood_ijk = log(binopdf(df_AB.g_ijk, df_AB.n_ijk, df_AB.p_hat_ij));
total_log_likelihood_AB = sum(df_AB.log_likelihood_ijk);

%% Row (A-only) model M_A: one probability per level of A
p_hat_A = varfun(@(x) sum(x), df, 'InputVariables', {'g_ijk', 'n_ijk'}, ...
    'GroupingVariables', {'Fac_A'});
p_hat_A.p_hat_ij = p_hat_A.Fun_g_ijk ./ p_hat_A.Fun_n_ijk;

df_A = innerjoin(df, p_hat_A(:, {'Fac_A', 'p_hat_ij'}), 'Keys', {'Fac_A'});
df_A.log_likelihood_ijk = log(binopdf(df_A.g_ijk, df_A.n_ijk, df_A.p_hat_ij));
total_log_likelihood_A = sum(df_A.log_likelihood_ijk);

%% Column (B-only) model M_B: one probability per level of B
p_hat_B = varfun(@(x) sum(x), df, 'InputVariables', {'g_ijk', 'n_ijk'}, ...
    'GroupingVariables', {'Fac_B'});
p_hat_B.p_hat_ij = p_hat_B.Fun_g_ijk ./ p_hat_B.Fun_n_ijk;

df_B = innerjoin(df, p_hat_B(:, {'Fac_B', 'p_hat_ij'}), 'Keys', {'Fac_B'});
df_B.log_likelihood_ijk = log(binopdf(df_B.g_ijk, df_B.n_ijk, df_B.p_hat_ij));
total_log_likelihood_B = sum(df_B.log_likelihood_ijk);

%% Additive model M_{A+B}: logit^{-1}(mu + alpha_i + beta_j), with identifiability via L
% The L matrix maps (mu, alpha_{2..I}, beta_{2..J}) to cell probabilities p_{ij}
% under additivity (no interaction). Here we follow the same construction as the R code.

% Build an I x I matrix with the first row as the reference and rows 2..I constrained
A_id = eye(I);
A_id(2:I) = 1;  % identify alpha with respect to level 1

% Construct L = [L1, L2]
L1 = kron(A_id,             ones(J,1));                 % A effects expansion
L2 = kron(ones(I,1), [zeros(1,J-1); eye(J-1)]);         % B effects expansion
L  = [L1, L2];                                          % maps theta -> cell logits (then to probabilities inside nlogApB_fun_un)

% Linear inequality constraints encode order/centering constraints used in the R code
constr_L = L(2:end, :);                                 % drop first row for baseline
Aineq = [constr_L; -constr_L];                          % A in fmincon (A*x <= b)
Bineq = [ones(size(constr_L, 1), 1); zeros(size(constr_L, 1), 1)]; % b in fmincon

% Initial values for theta = [mu, alpha_{2..I}, beta_{2..J}]'
% You can choose either an M_B-based or M_A-based initializer.
temp2_B = p_hat_B.p_hat_ij;
theta_ini_B = [temp2_B(1); zeros(I-1, 1); temp2_B(2:J) - temp2_B(1)];   % model-B based initial

temp2_A = p_hat_A.p_hat_ij;
theta_ini_A = [temp2_A(1); temp2_A(2:I) - temp2_A(1); zeros(J-1, 1)];   % model-A based initial

theta_ini = theta_ini_B;   % <- choose one; swap to theta_ini_A if preferred

% Parameter bounds (first element typically acts like an intercept on the link scale)
lb = [0; -Inf(I + J - 2, 1)];
ub = [1;  Inf(I + J - 2, 1)];

% Objective: negative log-likelihood under the additive model with unbalanced replicates
eval_f0 = @(theta) nlogApB_fun_un(theta, df, L, rep_vec);

% Fit M_{A+B} via constrained optimization
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display','off');
theta = fmincon(eval_f0, theta_ini, Aineq, Bineq, [], [], lb, ub, [], options);

% Predicted probabilities per replicate (expand each cell's p_{ij} by its replicate count)
p_hat_ij = repelem(L * theta, rep_vec);

df_ApB = df;
df_ApB.p_hat_ij = p_hat_ij;
df_ApB.log_likelihood_ijk = log(binopdf(df_ApB.g_ijk, df_ApB.n_ijk, p_hat_ij));
total_log_likelihood_ApB = sum(df_ApB.log_likelihood_ijk);

%% Likelihood ratio statistics (real-data, single run)
% Interaction: compare M_{AB} (saturated) vs M_{A+B} (additive)
log_diff_ApBvsAB = 2 * (total_log_likelihood_AB - total_log_likelihood_ApB);

% Factor A: compare M_{A+B} vs M_B
log_diff_BvsApB  = 2 * (total_log_likelihood_ApB - total_log_likelihood_B);

% Factor B: compare M_{A+B} vs M_A
log_diff_AvsApB  = 2 * (total_log_likelihood_ApB - total_log_likelihood_A);

%% P-values from the asymptotic chi-square reference
% df for interaction = (I-1)*(J-1)
% df for A main      = I-1
% df for B main      = J-1
Pval_interaction = chi2cdf(log_diff_ApBvsAB, (I-1)*(J-1), 'upper');  % H0: no interaction
Pval_A = chi2cdf(log_diff_BvsApB, (I-1), 'upper');                   % H0: no A effect
Pval_B = chi2cdf(log_diff_AvsApB, (J-1), 'upper');                   % H0: no B effect

% Display results
fprintf('Binomial LRT p-values:\n');
fprintf('  Interaction (A¡¿B): %g\n', Pval_interaction);
fprintf('  Factor A        : %g\n', Pval_A);
fprintf('  Factor B        : %g\n', Pval_B);







