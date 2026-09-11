function ss = fB0(n, A, I, J, adf, ab)
    IJ = I*J; B = zeros(IJ,1); t = 1;
    eps_ridge = 1e-12;
    mu = ab(1); alpha = ab(2:1+I);
    for i = 1:I
        for j = 1:J
            resid = A(1:n(t), t) - mu - alpha(i);
            sse = sum(resid.^2, 'omitnan');
            if ~isfinite(sse) || sse <= 0, sse = eps_ridge; end
            B(t) = - adf(t) * ((n(t)-1)/2) * log(sse + eps_ridge);
            t = t + 1;
        end
    end
    ss = sum(B);
end
