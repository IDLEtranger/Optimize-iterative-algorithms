function [value, gradient_q, hessian_q] = rosenbrock(x)
% input: [x; y] output: value gradient hessian
value = 100*( x(2) - x(1)^2 )^2 + (1 - x(1))^2;

if nargout >= 2 % need return gradient
    gradient_q = [(-400*x(1)*( x(2) - x(1)^2 ) - 2*( 1 - x(1)));
                200*( x(2) - x(1)^2 )];
    
    if nargout >= 3 % need return hessian matrix
        hessian_q = [(1200*x(1) - 400*x(2) +2), -400*x(1);
                   -400*x(1), 200];
    end

end