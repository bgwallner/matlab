clear 
syms k q2 m2 g alpha mu S q1 F

% Kraftj�mvikt sf�ren
%----------------------
% Friktionskraft motverkar r�relsen och sf�ren "r�r sig" 
% upp f�r kilen, dvs friktionskraft fr�n kil riktad snett 
% ned�t l�ngs med kilens ovansida. Friktionskraft fr�n 
% v�ggen riktad ned�t.
% Fj�derkraft F = k*q2.

% N1*cos(alpha)-mu*N1*sin(alpha)-m2*g-mu*N2-k*q2 = 0
% mu*N1*cos(alpha)+N1*sin(alpha)-N2 = 0

% N1 = (m2*g+k*q2)/(cos(alpha)*(1-mu^2))
% N2 = (m2*g+k*q2)*(sin(alpha)+mu*cos(alpha))/(cos(alpha)*(1-mu^2))

% Kraftj�mvikt kilen
%---------------------
% Friktionskraft fr�n sf�ren riktad snett upp�t l�ngs med kilens
% ovansida, dvs motverkar r�relsen orsakad av kraften S.

% S - N1*(sin(alpha)+mu*cos(alpha)) = 0

Eq = S - (F + m2*g)*(sin(alpha)+mu*cos(alpha))/(cos(alpha)*(1-mu^2));
Eq2 = solve(Eq == 0, F);
pretty(Eq2)

% Antag att m2 = 0.
m2 = 0;
Eq3 = subs(Eq2);
pretty(Eq3)

% Evaluera Eq3
mu = 0:0.01:1;
alpha = pi/3;
S = 1;
g = matlabFunction(Eq3)
y = feval(g,1,alpha,mu);
plot(mu,y)
grid on
hold on

alpha = pi/4;
y = feval(g,1,alpha,mu);
plot(mu,y)
grid on
hold on

alpha = pi/8;
y = feval(g,1,alpha,mu);
plot(mu,y)
grid on
hold on

legend('\pi/3','\pi/4','\pi/8')
title('Relation mellan fjaderkraft och applicerad kraft F = k*S')
xlabel('Friktionskoefficient \mu')
ylabel('k')


