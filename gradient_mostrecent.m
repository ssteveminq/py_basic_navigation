% %Parameters
clc 
clear all 
close all
K=50; %attractive potential 
E=8; %repulsive potential
% area_width=12; %potential area width (m)
sx=2; % start pos x
sy=2; %start pos y
gx=3; %goal pos x
gy=-4; %goal pos y
ox=[5 6]; %obstacle x pos
oy= [-3 8]; %obstacle y pos
% [x,y]= meshgrid(-12:.5:12); 
[x,y] = meshgrid(-12:1.0:12,-12:1.0:12);
r=((x-gx).^2+(y-gx).^2).^.5;%Distance to Goal
r1=((x-ox(1)).^2+(y-oy(1)).^2).^.5;%Distance to Obs1
r2=((x-ox(2)).^2+(y-oy(2)).^2).^.5;%Distance to Obs2
U=K.*r+E./r1*+E./r2;
figure(1)
%hold on
figure(1)
mesh(x,y,U);
figure(2)
contour(x,y,U);
[Ux,Uy]=gradient(U,1.0);

% [J]=-gradient(U);
% gradx=zeros(49,49);
% grady=zeros(49,49);
% i=0;
% j=0;
% for x0=-24:1:24
%     i=i+1;
%     for y0=-24:1:24
%         j=j+1;
%         t = (x == x0*.5) & (y == y0*.5);
%         indt = find(t);
%         U_grad = [Ux(indt) Uy(indt)];
%         gradx(i,j)=U_grad(1);
%         grady(i,j)=U_grad(2);
%       
%     end
% 
% end
% gradx=grad(:,1);
% grady=grad(:,2);
% hold on
[u,v,w] = surfnorm(x,y,U);
figure
% [x1,y1]=meshgrid(
% quiver(x,y,0.1*Ux,0.1*Uy)
quiver3(x,y,U,u,v,w,0.5)
% set(h1,'AutoScale','on', 'AutoScaleFactor', 0.5)
%plot(Ux,Uy);s
%quiver(x,y,Ux,Uy)
%plot(Ux,Uy);
xnew=0;
ynew=0;
i=1;
iter=0;
% x=-2;
% y=2;
gvecx=[];
gvecy=[];
% maxiter=100;
errx=1;
erry=1;
alpha=.01;
tol=1e-4;
xi=sx;
yi=sy;

%while errx>tol && erry>tol
for i=1:49 %used to say -12:.5:12 
   
      for j=1:49
            

            xnew=xi-alpha*Ux(i,j);
%             errx = abs((xnew-xi)/xi);
            ynew=yi-alpha*Uy(i,j);
%             err = abs((ynew-yi)/yi);

            gvecx=[gvecx;xi];
            gvecy=[gvecy;yi];

            iter=iter+1;
            j=j+1;
            xi=xnew;
            yi=xnew;
            
%         disp("x is")
%         disp(x)
%         disp("y is")
%         disp(y)
          
%         end
      end
     i=i+1;
end

figure(3)
%r=zeros(length(gvecx));
%plot3(gvecx,gvecy,r)
hold on 
scatter(gvecx,gvecy);
contour(x,y,U);
hold off
%hold off
% figure(2)
% hold on
% contour(x,y,U)
% quiver(x,y,Ux,Uy)
% hold off
