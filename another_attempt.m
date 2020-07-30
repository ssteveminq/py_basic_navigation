clc
clear all 
close all
alpha=.05;
Iter=1;
Maxiter=101;
xs=-4;
ys=-4;
gx=-1;
gy=-2;
xo1=5;
yo1=-3;
xo2=6;
yo2=6;
E=80;
K=50;
% U=@(xs,ys)((K*((xs-gx)^2+(ys-gy)^2)^.5)+(E/((xs-xo1)^2+(ys-yo1)^2)^.5)+(E/((xs-xo2)^2+(ys-yo2)^2)^.5));
% Udx=@(xs,ys)((diff(U(xs,ys),xs)));
% Udy=@(xs,ys)((diff(U(xs,ys),ys)));
Udx=@(xs,ys)((50*(xs-1))/(xs^2-2*xs+ys^2+5-4*ys)^.5)-((80*(xs-5))/(xs^2-10*xs+ys^2+34-6*ys)^1.5)-((80*(xs-6))/(xs^2-12*xs+ys^2+72-12*ys)^1.5);
Udy=@(xs,ys)((50*(ys-2))/(ys^2-4*ys+xs^2+5-2*xs)^.5)-((80*(ys-3))/(ys^2-6*ys+xs^2+34-10*xs)^1.5)-((80*(ys-6))/(ys^2-12*ys+xs^2+72-12*xs)^1.5);
while Iter<Maxiter && xs~=gx && ys~=gy
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