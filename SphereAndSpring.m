clear 
syms k q2 m2 g alpha mu S q1 F

% Stationary part of Lagranges equation
N2 = (k*q2 + m2*g)*(tan(alpha) - mu)/(1 + mu^2);
Fgen = -2*mu*N2*tan(alpha);
LH = m2*g*tan(alpha) + k*(q2/tan(alpha))*(tan(alpha))^2;

% Solve equation
Seq = solve(LH == S - Fgen, S);
pretty(Seq);

m2 = 0;
q2 = q1*tan(alpha);
Seq2 = subs(Seq);
pretty(Seq2);

q1 = F/(k*tan(alpha));
Seq3 = subs(Seq2);
pretty(Seq3);

% Solve equation
Seq4 = solve(S == Seq3, F);
pretty(Seq4)