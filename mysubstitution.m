clear all
clc

f = @(t) 3.402*10^(-9)*(t^4-300^4)
t0 = input('first value : ')
t1 = input('second value : ')
f0 = f(t0);
f1 = f(t1);
t2 = f0;
f2 = f(t2)
err = 0.001;

while f2 > err
    t0 = f2;
    f0 = f(t0);
    t2 = f0;
    f2 = f(t2);
end
disp (['temp :', num2str(t2)])

