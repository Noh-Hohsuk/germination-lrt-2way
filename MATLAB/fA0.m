function ss = fA0(n, A, I, J, adf, ab)
    % H0(A): mu + beta_j (A È¿°ú ¾øÀ½), no interaction
    % ab = [mu, beta(1:J)]
    IJ = I*J; B = zeros(IJ,1); t = 1;
    eps_ridge = 1e-12;
    mu = ab(1); beta = ab(2:1+J);
    for i = 1:I
        for j = 1:J
            resid = A(1:n(t), t) - mu - beta(j);
            sse = sum(resid.^2, 'omitnan');
            if ~isfinite(sse) || sse <= 0, sse = eps_ridge; end
            B(t) = - adf(t) * ((n(t)-1)/2) * log(sse + eps_ridge);
            t = t + 1;
        end
    end
    ss = sum(B);
end
