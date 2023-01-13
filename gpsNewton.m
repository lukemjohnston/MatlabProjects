function [v] = gpsNewton(sat1, sat2, sat3, sat4, x0, y0, z0, d0)

cc=299792.458;  %speed of light

A = [sat1(1), sat2(1), sat3(1), sat4(1)];
B = [sat1(2), sat2(2), sat3(2), sat4(2)];
C = [sat1(3), sat2(3), sat3(3), sat4(3)];
t = [sat1(4), sat2(4), sat3(4), sat4(4)];

v = [x0, y0, z0, d0];

for k=1:4
    x = v(1);
    y = v(2);
    z = v(3);
    d = v(4);

    F = [((x - A(1)).^2 + (y - B(1)).^2 + (z - C(1)).^2 - (cc*(t(1) - d)).^2);
    ((x - A(2)).^2 + (y - B(2)).^2 + (z - C(2)).^2 - (cc*(t(2) - d)).^2);
    ((x - A(3)).^2 + (y - B(3)).^2 + (z - C(3)).^2 - (cc*(t(3) - d)).^2);
    ((x - A(4)).^2 + (y - B(4)).^2 + (z - C(4)).^2 - (cc*(t(4) - d)).^2); ];

    J = [(2*(x-A(1))) (2*(y-B(1))) (2*(z-C(1))) (2*cc*(cc*t(1)-cc*d));
        (2*(x-A(2))) (2*(y-B(2))) (2*(z-C(2))) (2*cc*(cc*t(2)-cc*d));
        (2*(x-A(3))) (2*(y-B(3))) (2*(z-C(3))) (2*cc*(cc*t(3)-cc*d));
        (2*(x-A(4))) (2*(y-B(4))) (2*(z-C(4))) (2*cc*(cc*t(4)-cc*d)); ];

    H = J\F;

    v(1) = v(1) - H(1);
    v(2) = v(2) - H(2);
    v(3) = v(3) - H(3);
    v(4) = v(4) - H(4);

end

end