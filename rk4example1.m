% Runge-Kutta Method of order 2
% general formula for Runge-Kutta Method of order 2
%k1 = h*f(tn,yn);
%k2 = h*f(t(n+1),y(n)+k1);
%y(n+1) = y(n) + (1/2)*(k1 + k2);
% example one: motion of a falling object with air resistance satisfying
% the equation dv/dt = 9.8 - 0.2*v, v(0) = 0
f = @(t,v)9.8-0.2*v;
t0 = 0;
v0 = 0 ;
h = 2 ;
tn = 50 ;
n = (tn - t0)/h;
t(1) = t0; v(1) = v0;
for i = 1:n
    t(i+1) = t0 + i*h;
    k1 = h*f(t(i),v(i));
    k2 = h*f(t(i+1),v(i)+k1);
    v(i+1) = v(i) + (1/2)*(k1 + k2);
    fprintf('v(%.2f) = %.4f\n', t(i+1),v(i+1))
end
timetaken = toc;
disp(timetaken)
figure;
hold on
plot(t, v);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity of the Ball Over Time');
grid on;



