function [nll] = nlogApB_fun(theta,df,L,rep_num)

    p_hat_ij = kron(L * theta, ones(rep_num, 1));
    % Compute the log-likelihood
    log_likelihood_ijk = log(max(binopdf(df.g_ijk, df.n_ijk, p_hat_ij), 1e-100));

    % Return the negative mean log-likelihood
    nll=-mean(log_likelihood_ijk);
end