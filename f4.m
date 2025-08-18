function ss = f4(n, A, I, J, adf, ab)
    % H1: interaction 포함 모델의 통합로그우도(합) 계산
    % 안정화: sse가 0이 되는 경우 log(0) 회피를 위해 작은 ridge 추가

    IJ = I * J;
    B  = zeros(IJ, 1);          % 미리 할당
    t  = 1;
    eps_ridge = 1e-12;          % 수치 안정화(필요시 1e-10 등으로 조정)

    for i = 1:I
        for j = 1:J
            % 잔차: Y - (mu + a_i + b_j + gamma_ij)
            resid = A(1:n(t), t) - ab(1) - ab(i+1) - ab(I+j+1) - ab(I+J+1+t);

            % SSE 계산 (NaN 안전)
            sse = sum(resid.^2, 'omitnan');

            % sse가 0 또는 음수/비정상일 때 보호
            if ~isfinite(sse) || sse <= 0
                sse = eps_ridge;
            end

            % 통합로그우도 항(부호 주의: 원래 코드 구조 유지)
            B(t) = - adf(t) * ((n(t)-1)/2) * log(sse + eps_ridge);

            t = t + 1;
        end
    end

    ss = sum(B);
end
