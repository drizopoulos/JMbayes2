data {
    int n;
    int p;
    matrix[n, p] b;
    matrix[n, p] mu_b;
    real<lower=0> scale_diag_D;
    real<lower=0> lkj_shape;
}

parameters {
    vector<lower = 0>[p] L_var_D;
    cholesky_factor_corr[p] L_corr_D;
}

model {
    matrix[p, p] L_D;
    L_D = diag_pre_multiply(L_var_D, L_corr_D);
    L_var_D ~ student_t(3, 0, scale_diag_D);
    L_corr_D ~ lkj_corr_cholesky(lkj_shape);
    for (i in 1:n) {
        b[i, ] ~ multi_normal_cholesky(mu_b[i, ], L_D);
    }
}

generated quantities {
    matrix[p, p] D;
    D = diag_pre_multiply(L_var_D, L_corr_D) * diag_pre_multiply(L_var_D, L_corr_D)';
}
