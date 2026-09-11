function ss = f3(n, A, I, J, adf, ab)
    t = 1; B = zeros(I*J,1);
    eps_ridge = 1e-12;              
    for i = 1:I
        for j = 1:J
            resid = A(1:n(t), t) - ab(1) - ab(i+1) - ab(I+j+1);
            sse = sum(resid.^2);
            B(t) = -adf(t) * ((n(t)-1)/2) * log(sse + eps_ridge);   
            t = t + 1;
        end
    end
    ss = sum(B);
end
