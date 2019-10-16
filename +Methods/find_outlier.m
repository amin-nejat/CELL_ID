function outlier = find_outlier(aligned, mu, sigma)
ll_all = 2*Methods.AutoId.pdist2_maha(aligned, mu, sigma);
p_all = chi2cdf(ll_all,6,'upper');
ll_pos = 2*Methods.AutoId.pdist2_maha(aligned(1:3), mu(1:3), sigma(1:3,1:3));
p_pos = chi2cdf(ll_pos,3,'upper');
ll_col = 2*Methods.AutoId.pdist2_maha(aligned(4:6), mu(4:6), sigma(4:6,4:6));
p_col = chi2cdf(ll_col,3,'upper');
for i=1:3
    icol = i+3;
    ll_col_single = 2*Methods.AutoId.pdist2_maha(aligned(icol), mu(icol), sigma(icol,icol));
    p_col_single(i) = chi2cdf(ll_col_single,1,'upper');
    ll_pos_single = 2*Methods.AutoId.pdist2_maha(aligned(i), mu(i), sigma(i,i));
    p_pos_single(i) = chi2cdf(ll_pos_single,1,'upper');
end

outlier.global_test = p_all;
outlier.position_test.global = p_pos;
outlier.position_test.individual = p_pos_single;
outlier.color_test.global = p_col;
outlier.color_test.individual = p_col_single;
end