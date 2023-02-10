clear all
clc

R = input('enter the range : ')
f = @(v) (v*v/9.8) - R
v0 = input('first guess : ')
v1 = input('second guess : ')
err = 0.01;

f0 = f(v0);
f1 = f(v1);

if f0*f1 > 0
    disp(['not satisfy'])
else
    v2 = v0 - ((v1-v0)*f0/(f1-f0));
    f2 = f(v2);
    while abs(f2) > err 
        if f2*f1<0
            v0 = v2;
            f0 = f2;
        else
            v1 = v2;
            f1 = f2;
        end
       v2 = v0 - ((v1-v0)*f0/(f1-f0));
       f2 = f(v2);
    end
    disp(['velocity : ', num2str(v2)])
end