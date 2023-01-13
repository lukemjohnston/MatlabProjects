format shortG;

sat1 = [15600 7540 20140 0.07074];
sat2 = [18760 2750 18610 0.07220];
sat3 = [17610 14630 13480 0.07690];
sat4 = [19170 610 18390 0.07242];

x0=0;
y0=0;
z0=6700;
d0=0;

v = gpsNewton(sat1, sat2, sat3, sat4, x0, y0, z0, d0);
disp(v)