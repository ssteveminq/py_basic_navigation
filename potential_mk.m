close all
clear all

syms X Y
initial = [2,-2];
goal = [6,-2];
obs = [0,1];
obs2 = [-5,8];
K_att = 1000;
gamma = 900;
x=6;
y=-2;

[X,Y] = meshgrid(-12:.05:12);
Uatt = K_att*((goal(1)-X).^2 +(goal(2)-Y).^2);
Urep = (gamma*1./((obs(1)-X).^2 +(obs(2)-Y).^2)-200);
Urep2 = (gamma*1./((obs2(1)-X).^2 +(obs2(2)-Y).^2)-200);

U = 1000.*(6-X).^2+(-2-Y).^2+900./((-X).^2+(1-Y).^2-200)+900./((-5-X).^2+(8-Y).^2-200);
Upx=2000.*X - 1800.*X./(X.^2 + (1 - Y).^2 - 200).^2 + 900.*(-2.*X - 10)./((8 - Y).^2 + (-X - 5).^2 - 200).^2 - 12000;
Upy=2.*Y + 900.*(2 - 2.*Y)./(X.^2 + (1 - Y).^2 - 200).^2 + 900.*(16 - 2.*Y)./((8 - Y).^2 + (-X - 5).^2 - 200).^2 + 4;
Utotal = Uatt+Urep+Urep2;
%[Upx,Upy]=gradient(U);

tol=.001;
% [Upx,Upy]=gradient(Utotal);
%quiver(X,Y,Upx,Upy);
% subplot(3,1,1)
% mesh(X,Y,Uatt)
% subplot(3,1,2)
% mesh(X,Y,Urep+Urep2)
% subplot(3,1,3)
% mesh(X,Y,Utotal)

while x~=6 && y~=-2
    for i=1:200
        Un=Upx(i).*Upx(i)+Upy(i).*Upy(i);
        Upx2=Upx./Un;
        Upy2=Upy./Un;
        X(i)=X(i-1)-alpha*Upx2;
        Y(i)=Y(i-1)-alpha*Upy2;
        disp(X)
        disp(Y)
    end
    i=i+.01;
end
%plot(X,Y)
