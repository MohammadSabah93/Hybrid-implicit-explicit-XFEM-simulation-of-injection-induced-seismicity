function [S] = slipfunction(params, Ss_trial, slipfix, thetafix, Sn, dt, damp)

Ks=params.Ks;
Kn=params.Kn;
V0=params.V0;
Dc=params.Dc;
a=params.a;
b=params.b;
uf0=params.uf0;
S=slipfix;
iter = 0;
Norm = Inf;
tolerance = 1e-8;
maxIter = 100;
while Norm > tolerance && iter < maxIter

iter = iter + 1;


V = (S-slipfix)/dt;

if abs(V)<V0
  V=V0; 
end

theta = (thetafix+dt)/(1+(dt*abs(V))/Dc); 
uf = uf0 + a*log(abs(V)/V0) + b*log(V0*theta/Dc);
r = abs(Ss_trial - Ks * sign(Ss_trial) * (S-slipfix)) - uf * abs(Sn) - damp * V;
duf = 1/dt * (a/V - b/theta*(Dc*dt*(theta+dt))/((Dc+dt*V)^2));

J = Ks - duf * Sn + damp/dt;
S = S + (r/J);
if (abs(norm(r))<=tolerance)
    break;
end
end