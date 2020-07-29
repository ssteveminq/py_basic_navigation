clc
clear all 
close all

alpha=.005;
Iter=1;
Maxiter=5;
xs=10;
ys=10;
Udx=@(xs,ys)((50*(xs-1))/(xs^2-2*xs+ys^2+5-4*ys)^.5)-((80*(xs-5))/(xs^2-10*xs+ys^2+34-6*ys)^1.5)-((80*(xs-6))/(xs^2-12*xs+ys^2+72-12*ys)^1.5);
Udy=@(xs,ys)((50*(ys-2))/(ys^2-4*ys+xs^2+5-2*xs)^.5)-((80*(ys-3))/(ys^2-6*ys+xs^2+34-10*xs)^1.5)-((80*(ys-6))/(ys^2-12*ys+xs^2+72-12*xs)^1.5);
while Iter<Maxiter
    xs=xs-alpha*Udx(xs,ys);
    ys=ys-alpha*Udx(xs,ys);
    disp("Iter=");
    disp(Iter);
    disp("x=")
    disp(xs);
    disp("y=");
    disp(ys);
    Iter=Iter+1;
end