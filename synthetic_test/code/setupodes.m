function fout = setupodes( t, H, lambda, parameters, pressurestruct)
%Set up ODEs for H = [h; dh_dh0; dh_dlambda_k]
% fout(1) = dh/dt
% fout(2) = d/dt (h_h0)
% fout(3...kmax+2) = d/dt (h_lambda_k)
% We will pass fout into ode45(...) to solve for each of its components.
% This function is called in solveforH.m

h = H(1);
h_h0 = H(2);
if (length(H) > 2)
    h_lambda = H(3:end);
    fout = zeros(parameters.n_lambda + 2, 1);
else
    fout = zeros(2,1);
end

[F, dFdh, dFdlambda] = sheetthickness(t, H, lambda, parameters, pressurestruct);

fout(1) = F; %dh/dt = F

fout(2) = dFdh*h_h0; %d/dt (h_h0) = dF/dh * h_h0

if (length(H) > 2)
    fout(3:end) = dFdlambda + ...
        dFdh*h_lambda; %d/dt(h_lambda_k) = dF/dlambda_k + (dF/dh)*h_lambda_k
end

end


