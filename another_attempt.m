clc 
clear all
close all
alpha=.01;
Iter=1;
Maxiter=626;
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
U=@(xs,ys)((K*((xs-gx)^2+(ys-gy)^2)^.5)+(E/((xs-xo1)^2+(ys-yo1)^2)^.5)+(E/((xs-xo2)^2+(ys-yo2)^2)^.5));
% Udx=@(xs,ys)((diff(U(xs,ys),xs)));
% Udy=@(xs,ys)((diff(U(xs,ys),ys)));
Udx=@(xs,ys)((50*(xs-1))/(xs^2-2*xs+ys^2+5-4*ys)^.5)-((80*(xs-5))/(xs^2-10*xs+ys^2+34-6*ys)^1.5)-((80*(xs-6))/(xs^2-12*xs+ys^2+72-12*ys)^1.5);
Udy=@(xs,ys)((50*(ys-2))/(ys^2-4*ys+xs^2+5-2*xs)^.5)-((80*(ys-3))/(ys^2-6*ys+xs^2+34-10*xs)^1.5)-((80*(ys-6))/(ys^2-12*ys+xs^2+72-12*xs)^1.5);
dx=zeros(25,25);
dy=zeros(25,25);
for i=1:25
    for j=1:25
        xs=-12+1*(i-1);
        ys=-12+1*(j-1);
        dx(i,j)=Udx(xs,ys);
        dy(i,j)=Udy(xs,ys);
%         xs=xs-alpha*Udx(xs,ys);
%         ys=ys-alpha*Udx(xs,ys);
   
        disp("Iter=");
        disp(Iter);
        disp("x=")
        disp(xs);
        disp("y=");
        disp(ys);
        Iter=Iter+1;
    end
end

% while xs~=gx && ys~=gy
%     dx=[dx;Udx(xs,ys)];
%     dy=[dy;Udy(xs,ys)];
%     xs=xs-alpha*Udx(xs,ys);
%     ys=ys-alpha*Udx(xs,ys);
%    
%     disp("Iter=");
%     disp(Iter);
%     disp("x=")
%     disp(xs);
%     disp("y=");
%     disp(ys);
%     Iter=Iter+1;
% end
%%
clc 
close all 
clear all
syms  x y
U=@(xs,ys)((K*((xs-gx)^2+(ys-gy)^2)^.5)+(E/((xs-xo1)^2+(ys-yo1)^2)^.5)+(E/((xs-xo2)^2+(ys-yo2)^2)^.5));

%%
clc
clear all
close all
K=50; %attractive potential 
E=80; %repulsive potential

sx=10; % start pos x
sy=10; %start pos y
gx=-1; %goal pos x
gy=-2; %goal pos y
ox=[5 6]; %obstacle x pos
oy= [-3 6]; %obstacle y pos 

[a,b] = meshgrid(-12:.5:12,-12:.5:12);
r=((a-gx).^2+(b-gx).^2).^.5;%Distance to Goal
r1=((a-ox(1)).^2+(b-oy(1)).^2).^.5;%Distance to Obs1
r2=((a-ox(2)).^2+(b-oy(2)).^2).^.5;%Distance to Obs2
Pot=K.*r+E./r1*+E./r2;
hold on
contour(a,b,Pot);


xs=10;
ys=10;
gx=-1;
gy=-2;
xo1=5;
yo1=-3;
xo2=6;
yo2=6;
E=80;
K=50;
alpha=.05;
% syms x y 
% U=@(x,y)((K*((x-gx)^2+(y-gy)^2)^.5)+(E/((x-xo1)^2+(y-yo1)^2)^.5)+(E/((x-xo2)^2+(y-yo2)^2)^.5));
% disp(pdx(1,1));
% disp(pdy(1,1));
xpos=[10];
ypos=[10];
iter=0;
maxiter=100;
tol=1e-1;
errx=5;
erry=5;
xsn=0;
ysn=0;
while iter<maxiter && errx>tol && erry>tol
    
    xsn=xs-alpha*pdx(xs,ys);
    errx = abs((xsn-gx)/gx);
    ysn=ys-alpha*pdy(xs,ys);
    erry= abs((ysn-gy)/gy);
    xpos=[xpos;xsn];
    ypos=[ypos;ysn];
    iter=iter+1;
    xs=xsn;
    ys=ysn;
end
scatter(gx,gy,'r')
scatter(sx,sy,'g')
line(xpos,ypos)

% plot(xpos,ypos);
hold off
function[Udx]= pdx(xs,ys)
gx=-1;
gy=-2;
xo1=5;
yo1=-3;
xo2=6;
yo2=6;
E=80;
K=50;
syms x y 
U=@(x,y)((K*((x-gx)^2+(y-gy)^2)^.5)+(E/((x-xo1)^2+(y-yo1)^2)^.5)+(E/((x-xo2)^2+(y-yo2)^2)^.5));
pdx=diff(U,x);
Ux=subs(pdx,{x,y},{xs,ys});
Udx=double(Ux);
end 

function[Udy]= pdy(xs,ys)
gx=-1;
gy=-2;
xo1=5;
yo1=-3;
xo2=6;
yo2=6;
E=80;
K=50;
syms x y 
U=@(x,y)((K*((x-gx)^2+(y-gy)^2)^.5)+(E/((x-xo1)^2+(y-yo1)^2)^.5)+(E/((x-xo2)^2+(y-yo2)^2)^.5));
pdy=diff(U,y);
Uy=subs(pdy,{x,y},{xs,ys});
Udy=double(Uy);
end 

