function plot_cost_function_wrt_2_vars( switchingevents, pressurestruct, ...
    Beta, parameters, index1, index2, range1, range2, meshsize)
% Plot cost function C with respect to param 1 = Beta(index1) and 
%  param 2 = Beta(index2).
%  Domain for param 1 is [min1 max1] and domain for param 2 is [min2 max2].

if (nargin <= 7) %range1 and range2 not specified
    range1 = [Beta(index1) - 10, Beta(index1) + 10];
    range2 = [Beta(index2) - 10, Beta(index2) + 10];
end
if (nargin == 6 || nargin == 8) %meshsize not specified
    meshsize = .1;
end

a1 = range1(1);
b1 = range1(2);
a2 = range2(1);
b2 = range2(2);

param1 = a1:meshsize:b1;
param2 = a2:meshsize:b2;
[X,Y] = meshgrid(param1, param2); % X is meshgrid for param1, Y is meshgrid for param2
[n,m] = size(X);
C = 0*Y;
for i = 1:n
    for j = 1:m
        Beta(index1) = X(i,j);
        Beta(index2) = Y(i,j);
        [t, H] = solveforH(Beta(parameters.n_lambda+1), parameters, ...
            Beta(1:parameters.n_lambda), pressurestruct);
            %[t,H] = solveforH(h0, parameters, lambda, pressurestruct)
        disp('Solved for H.');
        h_SE = interp1(t, H(:,1), switchingevents, 'spline');
        disp('Solved for h_SE.');
        C(i,j) = calculatecostfunction(Beta, h_SE, parameters);
        disp(['Solved for C(' num2str(i) ',' num2str(j) ') = ' num2str(C(i,j))]);
    end
end


% if (index1 == parameters.h_r.index)
%     h_r = a1:meshsize:b1; % h_r is first independent variable
%     label1 = 'h_r';
% elseif (index2 == parameters.h_r.index)
%     h_r = a2:meshsize:b2; % h_r is second independant variable
%     label2 = 'h_r';
% else
%     h_r = Beta(parameters.h_r.index); % treat h_r as constant
% end
% 
% if (index1 == parameters.logk.index)
%     logk = a1:meshsize:b1; % logk is first independent variable
%     label1 = 'logk';
% elseif (index2 == parameters.logk.index)
%     logk = a2:meshsize:b2; % logk is second independant variable
%     label2 = 'logk';
% else
%     logk = Beta(parameters.logk.index); % treat logk as constant
% end
% 
% if (index1 == parameters.p_ice.index)
%     p_ice = a1:meshsize:b1; % p_ice is first independent variable
%     label1 = 'p_ice';
% elseif (index2 == parameters.p_ice.index)
%     p_ice = a2:meshsize:b2; % p_ice is second independant variable
%     label2 = 'p_ice';
% else
%     p_ice = Beta(parameters.p_ice.index); % treat p_ice as constant
% end
% 
% if (index1 == parameters.n_G.index)
%     n_G = a1:meshsize:b1; % n_G is first independent variable
%     label1 = 'n_G';
% elseif (index2 == parameters.n_G.index)
%     n_G = a2:meshsize:b2; % n_G is second independant variable
%     label2 = 'n_G';
% else
%     n_G = Beta(parameters.n_G.index); % treat n_G as constant
% end
% 
% if (index1 == parameters.n_lambda+1)
%     h0 = a1:meshsize:b1;
%     label1 = 'h_0';
% elseif (index2 == parameters.n_lambda+1)
%     h0 = a2:meshsize:b2;
%     label2 = 'h_0';
% else
%     h0 = Beta(parameters.n_lambda+1);
% end
% 
% if (index1 == parameters.n_lambda+2)
%     hc = a1:meshsize:b1;
%     label1 = 'h_c';
% elseif (index2 == parameters.n_lambda+2)
%     hc = a2:meshsize:b2;
%     label2 = 'h_c';
% else
%     hc = Beta(parameters.n_lambda+2);
% end

mesh(C);


end

