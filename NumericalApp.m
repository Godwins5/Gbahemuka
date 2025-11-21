%Differential Problem: Newton-Raphson
f1 = @(x) x.^3 - x - 2;
diffMethod = DifferentialMethods;
diffMethod.f = f1;
diffMethod.tol = 1e-4;
diffMethod.solve();

%Integral Problem: Runge-Kutta 2nd Order 
f2 = @(t, v) 9.8 - 0.2*v;
intMethod = IntegralMethods;
intMethod.f = f2;
intMethod.tol = 1e-4;
intMethod.solve();