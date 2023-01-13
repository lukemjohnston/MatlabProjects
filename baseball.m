
% First Pitch
x = -2.509;
y = 50;
z = 5.928;

Vx = 9.182;
Vy = -132.785;
Vz = -10.967;

Ax = -19.268;
Ay = 30.713;
Az = -16.580;

K = 0.005152949;

Cd = 0.392648;
Cl = 0.255819;

phi = 236.0038*pi/180;
g = 32.174;


% Second Pitch
% x = -2.43;
% y = 50;
% z = 6.46;
% 
% Vx = 9.46;
% Vy = -143.17;
% Vz = -9.5;
% 
% Ax = -23.08;
% Ay = 34.2;
% Az = -26.09;
% 
% K = 0.005316103;
% 
% Cd = 0.392648;
% Cl = 0.255819;
% 
% phi = 4.591151161;
% g = 32.174;


getTracjectory(x, y, z, Vx, Vy,Vz, Ax, Ay, Az)  
rkAN(x, y, z, Vx, Vy, Vz, K, Cd, Cl, phi, g)

function getTracjectory(x, y, z, Vx, Vy,Vz, Ax, Ay, Az)   
    n = 50;
    t = 0;
    xM = zeros(1,n);
    yM = zeros(1,n);
    zM = zeros(1,n);
    
    % Get Position
    for i=1:n
        xM(i) = x + (t*Vx) + (1/2)*(t^2)*Ax;
        yM(i) = y + (t*Vy) + (1/2)*(t^2)*Ay;
        zM(i) = z + (t*Vz) + (1/2)*(t^2)*Az;
        t = t+0.01;
        %fprintf("%f)\t%f \t %f \t %f\n", t, xM(i), yM(i), zM(i))
    end

    % Get Velocity
    accuray = 0.001;
    n2 = 50/accuray;

    vM = zeros(1, n2);
    tM = zeros(1, n2);

    i = 1;
    for y = 0 : 0.001 : 50
        [vM(i), tM(i)] = velocity(y, Vx, Vy,Vz, Ax, Ay, Az);
        i = i+1;
    end

    fprintf("Sportvison Method\nfinal positon values:\n      X         Y          Z\n")
    fprintf("%f\t%f \t %f \t %f \n", xM(n), yM(n), zM(n))


    tiledlayout(2,1)

    nexttile
    plot(tM, vM, 'b');
    xlabel('Time (Sec)');
    ylabel('Velocity ft/sec');
    title("Velocity")

    nexttile
    plot3(xM,yM,zM)
    title("Position")
end

function[v, t] = velocity(yAt, Vx, Vy, Vz, Ax, Ay, Az)
    % Get Velocity
    y0 = 50;
    t = (-Vy-sqrt(Vy^2-2*Ay*(y0-yAt)))/Ay;
    
    x = Vx + Ax*(t);
    z = Vz + Az*(t);
    y = Vy + Ay*(t);

    v = sqrt(x^2+z^2+y^2);
end


function rkAN(x, y, z, Vx, Vy,Vz, K, Cd, Cl, phi, g)

    % Setup Equations for ODE - This is the problem spot
    xEq = @(Vx, Vy, Vz) -K*Cd*(sqrt(Vx^2+Vy^2+Vz^2))*Vx - K*Cl*(sqrt(Vx^2+Vy^2+Vz^2))*Vy*sin(phi);
    yEq = @(Vx, Vy, Vz) -K*Cd*(sqrt(Vx^2+Vy^2+Vz^2))*Vy + K*Cl*(sqrt(Vx^2+Vy^2+Vz^2))*(Vx*sin(phi)-Vz*cos(phi));
    zEq = @(Vx, Vy, Vz) -K*Cd*(sqrt(Vx^2+Vy^2+Vz^2))*Vz + K*Cl*(sqrt(Vx^2+Vy^2+Vz^2))*Vy*cos(phi) - g;
    

    n = 100000;
    h = (0.5)/n;

    % Setup Matrices
    xM = zeros(1, n);
    yM = zeros(1, n);
    zM = zeros(1, n);

    vM = zeros(1, n);
    vXM = zeros(1, n);
    vYM = zeros(1, n);
    vZM = zeros(1, n);

    tM = zeros(1, n);

    % Set Intial Values
    xM(1) = x;
    yM(1) = y;
    zM(1) = z;

    vM(1) = sqrt(Vx^2 + Vy^2 + Vz^2);
    vXM(1) = Vx;
    vYM(1) = Vy;
    vZM(1) = Vz;

    tM(1) = 0;


    fprintf("\n\nNathan's Method\n")
    %fprintf("\n\nPitch Position Data\n")
    %fprintf("%f)\t%f \t %f \t %f \t \t V=%f \n", 0, vXM(1), vYM(1), vZM(1), vM(1))

    for i = 1:n
        % For X
        k1X = h * xEq(vXM(i), vYM(i), vZM(i));
        k2X = h * xEq(vXM(i) + k1X/2, vYM(i) + k1X/2, vZM(i) + k1X/2);
        k3X = h * xEq(vXM(i) + k2X/2, vYM(i) + k2X/2, vZM(i) + k2X/2);
        k4X = h * xEq(vXM(i) + k3X, vYM(i) + k3X, vZM(i) + k3X);
        vXM(i+1) = vXM(i) + (k1X + 2*k2X + 2*k3X + k4X) / 6;

        % For Y
        k1Y = h * yEq(vXM(i), vYM(i), vZM(i));
        k2Y = h * yEq(vYM(i) + k1Y/2, vYM(i)+ k1Y/2, vZM(i)+ k1Y/2);
        k3Y = h * yEq(vXM(i) + k2Y/2, vYM(i) + k2Y/2, vZM(i) + k2Y/2);
        k4Y = h * yEq(vXM(i) + k3Y, vYM(i) + k3Y, vZM(i) + k3Y);
        vYM(i+1) = vYM(i) + (k1Y + 2*k2Y + 2*k3Y + k4Y) / 6;

        % For Z
        k1Z = h * zEq(vXM(i), vYM(i), vZM(i));
        k2Z = h * zEq(vXM(i) + k1Z/2, vYM(i)+ k1Z/2, vZM(i)+ k1Z/2);
        k3Z = h * zEq(vXM(i) + k2Z/2, vYM(i) + k2Z/2, vZM(i) + k2Z/2);
        k4Z = h * zEq(vXM(i) + k3Z, vYM(i) + k3Z, vZM(i) + k3Z);
        vZM(i+1) = vZM(i) + (k1Z + 2*k2Z + 2*k3Z + k4Z) / 6;


        % Update Position
        xM(i+1) = xM(i) + vXM(i) * h;
        yM(i+1) = yM(i) + vYM(i) * h;
        zM(i+1) = zM(i) + vZM(i) * h;

        % Keep Track of Time and Velocity
        tM(i+1) = tM(i)+h;
        vM(i+1) = sqrt(vXM(i+1)^2+vYM(i+1)^2+vZM(i+1)^2);

        %fprintf("%f)\t%f \t %f \t %f \t \t V=%f \n", tM(i+1), xM(i+1), yM(i+1), zM(i+1), vM(i+1))
    end
    
    fprintf("final positon values:\n      X         Y           Z           V\n")
    fprintf("%f\t%f \t %f \t %f \n", xM(n+1), yM(n+1), zM(n+1), vM(n+1))


    tiledlayout(2,1)

    nexttile
    plot(tM, vM, 'b');
    xlabel('Time (Sec)');
    ylabel('Velocity ft/sec');
    title("Nathan's Method: Velocity")
    
    nexttile
    plot3(xM,yM,zM)
    title("Nathan's Method: Position")
end
