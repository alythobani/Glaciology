% F = u_t.*h_r - exp(logk).*h.*(abs(p_ice - p_w_t)).^(n_G-1).*(p_ice - p_w_t);

u_t = .1; % value for velocity used in sheetthickness.m
h_rval = (6.1344e-19)*(4e5)^3/0.1; % 
logkval = log(6.1344e-19);
h = 70000;
p_iceval = 916.7*9.80665*83.8;
p_w_t = 900*800;
n_Gval = 3;


Func_h_r = @(h_r) u_t.*h_r - exp(logkval).*h.*(abs(p_iceval - p_w_t)).^(n_Gval-1).*(p_iceval - p_w_t);
Func_logk = @(logk) u_t.*h_rval - exp(logk).*h.*(abs(p_iceval - p_w_t)).^(n_Gval-1).*(p_iceval - p_w_t);


