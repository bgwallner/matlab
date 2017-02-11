% Bo-Göran Wallner 2016.
% Calculating t = t(h) for draining of waterfilled tube with 
% an overpressure p. 
% Compare with https://mathtab.com/app_id=4136 where calculation 
% can be made without any overpressure.

clear 
syms t a2 h a1 theta1 theta2 td Pi Rho Atot1 Atot2 V1 V2 L Cl od g p pback R

% Area of circular segment (see e.g. Wolfram Mathworld)
A = R^2*acos((R-h)/R)-(R-h)*sqrt(2*R*h-h^2); 
V = L*A; m = Rho*V;

% dV/dt = dV/dh * dh/dt = flow by physics 
% dt = dV/dh * (1/flow) * dh 

% Integration gives t and many integrations give t = t(h)
% Calculate dV/dh 
dvdh = diff(V,h);

% Calculate flow * area with respect to h % pback static here but will increase with time!!!!! 
den = Cl*Pi*od^2/4*sqrt(2*g*h + 2*(p-pback)/Rho);

% Show symbolic expressions to be used if not % having a symbolic toolbox. 
integrand = -simplify(dvdh/den)

% From now on we do calculations numeric 
td = 2; R = td/2; 
od = 0.1; 
L = 5; 
p = 10^5; % 1 bar overpressure 
g = 9.81; 
Rho = 1000; 
Cl = 0.6; 
Pi = 3.1415; 
pback = 0; %128000; % 1 bar pressure acting back through orifice

% Evaluate with values 
eq1 = -subs(simplify(dvdh/den));
% Convert eq1 and eq2 to Matlab functions instead 
% of symbolic and create function handles. 
fcn1 = matlabFunction(eq1);

% Calculate integral (small imaginary parts due 
% to numerical problem with precision. 
steps = 200; 
for i=0:1:steps    
    upper = td - i*td/steps;
    time(i+1) = real(integral(fcn1,upper,0));    
    height(i+1)= upper;    
    volume(i+1)=L*R^2*acos((R-upper)/R)-(R-upper)*sqrt(2*R*upper-upper^2); 
end

% Plot results. 
subplot(2,1,1); 
plot(time, fliplr(height)); 
grid on 
hold on 
title('Water level') 
xlabel('Time [seconds]') 
ylabel('Height of water [meter]')

subplot(2,1,2); 
plot(time, fliplr(volume)); 
grid on 
title('The remaining volume') 
xlabel('Time [seconds]') 
ylabel('Volume [m^3]')