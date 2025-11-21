%% % Runge-Kutta Method of order 2
% general formula for Runge-Kutta Method of order 2
%k1 = h*f(tn,yn);
%k2 = h*f(t(n+1),y(n)+k1);
%y(n+1) = y(n) + (1/2)*(k1 + k2);
% Example two: cooling of a hot metal rod, following the equation dT/dt =
% -0.07(T-25)
f = @(t,T) -0.07*(T-25);
t0 = 0;
T0 = 200 ;
h = 10 ;
tn = 150 ;
n = (tn - t0)/h;
t(1) = t0; T(1) = T0;
for i = 1:n
    t(i+1) = t0 + i*h;
    k1 = h*f(t(i),T(i));
    k2 = h*f(t(i+1),T(i)+k1);
    T(i+1) = T(i) + (1/2)*(k1 + k2);
    fprintf('T(%.2f) = %.4f\n', t(i+1),T(i+1))
end
timetaken1 = toc;
disp(timetaken1)
figure;
plot(t, T);
xlabel('Time (s)');
ylabel('Temperature (0C)');
title('cooling of a hot metal rod');
grid on;
