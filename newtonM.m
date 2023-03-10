function [x] = newtonM(X)
x = X.';

for k=1:10
    F = [(x(1) + x(2) + x(3) -3);
        ((x(1))^2 + (x(2))^2 + (x(3))^2 -5);
        (exp(x(1)) + (x(1)*x(2)) - (x(1)*x(3)) -1)];

    J = [ 1 1 1;
        (2*x(1)) (2*x(2)) (2*x(3));
        (exp(x(1)) + x(2) - x(3)) x(1) (-x(1))];

    H = J\F;
    x = x-H;

end

end
