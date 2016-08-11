function [ dFuncdlambdak ] = calculate_divided_difference_approximation( Func, lambdak, delta )
% Approximate dFunc/dlambdak using divided difference approximation:
%
% dFunc/dlambdak approx. = [Func(lambdak + delta) - Func(lambdak)] / delta
%
% This should return a result that is accurate to within an error of order
%   sup(d^2Func/dlambdak^2 * delta/2)

dFuncdlambdak = (Func(lambdak+delta) - Func(lambdak))./delta;



end

