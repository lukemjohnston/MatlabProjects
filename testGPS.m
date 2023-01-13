format shortG;

%satellite locations and time corrections
sat1 = [9261 1214 24874 0.07407];
sat2 = [12820 -9599 21201 0.07055];
sat3 = [-10317 9597 22526 0.07789];
sat4 = [530 -2397 1120 0.07100];

%initial values
x0=0;
y0=0;
z0=15000;
d0=0;

v = gpsNewton(sat1, sat2, sat3, sat4, x0, y0, z0, d0);
disp(v)

%find the latitude and longitude
r = sqrt((v(1))^2+(v(2))^2+(v(3))^2);
long = rad2deg(atan2(v(2),v(1)));
lat = rad2deg(pi/2-acos(v(3)/r));
disp([lat, long]);
