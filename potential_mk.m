close all
clear all

syms X Y
initial = [2,-2];
goal = [6,-2];
obs = [0,1];
obs2 = [-5,8];
K_att = 1000;
gamma = 900;


[X,Y] = meshgrid(-12:.05:12);
Uatt = K_att*((goal(1)-X).^2 +(goal(2)-Y).^2);
Urep = (gamma*1./((obs(1)-X).^2 +(obs(2)-Y).^2)-200);
Urep2 = (gamma*1./((obs2(1)-X).^2 +(obs2(2)-Y).^2)-200);

U(X,Y)=@(X,Y)(1000.*(6-X).^2+(-2-Y).^2+900./((-X).^2+(1-Y).^2-200)+900./((-5-X).^2+(8-Y).^2-200));

Utotal = Uatt+Urep+Urep2;
[Upx,Upy]=gradient(U);
% [Upx,Upy]=gradient(Utotal);
% quiver(X,Y,Upx,Upy);
% 
% subplot(3,1,1)
% mesh(X,Y,Uatt)
% subplot(3,1,2)
% mesh(X,Y,Urep+Urep2)
% subplot(3,1,3)
% mesh(X,Y,Utotal)

