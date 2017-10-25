% Auralius Manurung : manurunga@yandex.com
%
% Collection of solvers with fixed steps
% For an ODE: y_dot = f(t,y)

clc
clear all
close all

% ------ Test 1 ------
[t, ya] = euler(@myode1, 1, 0, 0.02, 0.001);
[t, yb] = rk4(@myode1, 1, 0, 0.02, 0.001);
[t, yc] = rk38(@myode1, 1, 0, 0.02, 0.001);
[t, yd] = heun(@myode1, 1, 0, 0.02, 0.001);
figure
hold on
plot(t, ya, 'b')
plot(t, yb, 'r')
plot(t, yc, 'm')
plot(t, yd, 'k')
legend('Euler', 'RK4', 'RK3/8', 'Heun')
title('Stiff ODE, $\dot{y} = -1000y$', 'interpreter', 'latex');

% ------ Test 2 ------
[t, ya] = euler(@myode2, [1 1 1]', 0, 20, 0.001);
[t, yb] = rk4(@myode2, [1 1 1]', 0, 20, 0.001);
[t, yc] = rk38(@myode2, [1 1 1]', 0, 20, 0.001);
[t, yd] = heun(@myode2, [1 1 1]', 0, 20, 0.001);
figure
hold on
plot3(ya(1,:), ya(2,:), ya(3,:), 'b');
plot3(yb(1,:), yb(2,:), yb(3,:), 'r');
plot3(yc(1,:), yc(2,:), yc(3,:), 'm');
plot3(yd(1,:), yd(2,:), yd(3,:), 'k');
legend('Euler', 'RK4', 'RK3/8', 'Heun')
title('Lorenz, $\sigma = 10, \beta = 8/3, \rho = 28$', 'interpreter', 'latex');

% ------ Test 3 ------
[t, ya] = euler(@myode3, [2 0]', 0, 1000, 0.001);
[t, yb] = rk4(@myode3, [2 0]', 0, 1000, 0.001);
[t, yc] = rk38(@myode3, [2 0]', 0, 1000, 0.001);
[t, yd] = heun(@myode3, [2 0]', 0, 1000, 0.001);
figure
hold on
plot(ya(1,:), ya(2,:), 'b');
plot(yb(1,:), yb(2,:), 'r');
plot(yc(1,:), yc(2,:), 'm');
plot(yd(1,:), yd(2,:), 'k');
legend('Euler', 'RK4', 'RK3/8', 'Heun')
title('Van der Pol, $\mu = 100$', 'interpreter', 'latex');

% ------ Test 3 ------
[t, ya] = euler(@myode4, -1, 0, 10, 0.001);
[t, yb] = rk4(@myode4, -1, 0, 10, 0.001);
[t, yc] = rk38(@myode4, -1, 0, 10, 0.001);
[t, yd] = heun(@myode4, -1, 0, 10, 0.001);
figure
hold on
plot(t, ya, 'b');
plot(t, yb, 'r');
plot(t, yc, 'm');
plot(t, yd, 'k');
legend('Euler', 'RK4', 'RK3/8', 'Heun')
title('$\dot{y} = y * sin(t)$', 'interpreter', 'latex');

% ------ ODE Example 1 ------
function ydot = myode1(t,y)
    ydot = -1000*y;
end

% ------ ODE Example 2 ------
function ydot = myode2(t,y)
    sigma = 10;
    beta = 8/3;
    rho = 28;
    ydot = [sigma * (y(2) - y(1)); 
            y(1) * (rho - y(3)) - y(2);
            y(1) * y(2) - beta * y(3)];
end

% ------ ODE Example 3 ------
function ydot = myode3(t,y)
    Mu = 100;
    ydot = [y(2); Mu*(1-y(1)^2)*y(2)-y(1)];
end

% ------ ODE Example 4 ------
function ydot = myode4(t,y)
    ydot = y*sin(t);
end


% -------------------------------------------------------------------------

function [t, y] = euler(odefun, y0, tstart, tfinal, dt)
    t = tstart:dt:tfinal;
    y = zeros(length(y0),length(t));
    y(:,1) = y0;
    for k = 2 : length(t)
        ydot = odefun(t(k), y(:,k-1));
        y(:,k) = y(:,k-1)+ydot.*dt;
    end
end

% -------------------------------------------------------------------------

function [t, y] = rk4(odefun, yn, tstart, tfinal, dt)
% Based on http://www.aip.de/groups/soe/local/numres/bookcpdf/c16-1.pdf
% See page 711 Equ. 16.1.3

    t = tstart:dt:tfinal;
    y = zeros(length(yn),length(t));
    y(:,1) = yn;
    
    for k = 2 : length(t)
        yn = y(:,k-1);
        tn = t(k-1);
        k1 = dt * odefun(tn, yn);   
        k2 = dt * odefun(tn + dt / 2 , yn + k1 / 2);
        k3 = dt * odefun(tn + dt / 2, yn + k2 / 2);
        k4 = dt * odefun(tn + dt, yn + k3);
        
        y(:,k) = yn + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6;
    end
end

% -------------------------------------------------------------------------

function [t, y] = rk38(odefun, yn, tstart, tfinal, dt)
    t = tstart:dt:tfinal;
    y = zeros(length(yn),length(t));
    y(:,1) = yn;
    
    for k = 2 : length(t)
        yn = y(:,k-1);
        tn = t(k-1);
        k1 = dt * odefun(tn, yn);   
        k2 = dt * odefun(tn + dt / 3, yn + k1 / 3);
        k3 = dt * odefun(tn + dt * 2 / 3, yn + - k1 / 3 + k2);
        k4 = dt * odefun(tn + dt, yn + k1 - k2 + k3);
        
        y(:,k) = yn + 1 / 8 * k1 + 3 / 8 * k2 + 3 / 8 * k3 + 1 / 8 * k4;
    end
end

% -------------------------------------------------------------------------

function [t, y] = heun(odefun, yn, tstart, tfinal, dt)
    t = tstart:dt:tfinal;
    y = zeros(length(yn),length(t));
    y(:,1) = yn;
    
    for k = 2 : length(t)
        yn = y(:,k-1);
        tn = t(k-1);
        k1 = dt * odefun(tn, yn);   
        k2 = dt * odefun(tn + dt, yn + k1);
        
        y(:,k) = yn + 1 / 2 * k1 + 1 / 2 * k2;
    end
end
