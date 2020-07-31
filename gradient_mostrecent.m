% %Parameters
clc 
clear all 
close all
K=35; %attractive potential 
E=15; %repulsive potential
% area_width=12; %potential area width (m)
sx=10; % start pos x
sy=11; %start pos y
gx=-1; %goal pos x
gy=-5; %goal pos y
ox=[5.0 6.05]; %obstacle x pos
oy= [-3.0 6.05]; %obstacle y pos
% [x,y]= meshgrid(-12:.5:12); 
res=0.25;
[x,y] = meshgrid(-12:res:12,-12:res:12);
r=((x-gx).^2+(y-gy).^2).^.5;%Distance to Goal
r1=((x-ox(1)).^2+(y-oy(1)).^2).^.5;%Distance to Obs1
r2=((x-ox(2)).^2+(y-oy(2)).^2).^.5;%Distance to Obs2
U=K.*r+E./r1*+E./r2;

%hold on
figure(1)
mesh(x,y,U);
figure(2)
contour(x,y,U);
[Ux,Uy]=gradient(U,res,res);
% Ux=diff(U,x,1.0);
% Uy=diff(U,y,1.0);
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
figure(3)
% [x1,y1]=meshgrid(
% quiver(x,y,0.1*Ux,0.1*Uy)
% quiver3(x,y,U,u,v,w,0.5)
quiver(x,y,u,v,res)

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

% maxiter=100;
errx=1;
erry=1;
alpha=0.03;
tol=0.3;
xi=sx;
yi=sy;
gvecx=[];
gvecy=[];
%while errx>tol && erry>tol
% while xi~=gx && yi~=gy
%%
figure(4)
hold on
contour(x,y,U)
scatter(sx,sy,'b')
scatter(gy,gy,'r')

xmin =-12.0;
ymin =-12.0;
maxiter=500;

% for i=1:25 %used to say -12:.5:12 
   
%       for j=1:25
Ux=Ux';
Uy=Uy';
  
while iter<maxiter         
          ind_x = floor((xi-xmin)/res)
          ind_y = floor((yi-ymin)/res)
%           xx = xmin+res*ind_x
%           yy = ymin+res*ind_y
          gradx = Ux(ind_x,ind_y)
          grady = Uy(ind_x,ind_y)
%           if abs(gradx)<1.0
%              alpha = 0.2;
%           else
%               alpha= 0.05;
%           end
%           
%            if abs(grady)<1.0
%              alpha_y = 0.5;
%           else
%               alpha_y= 0.05;
%           end
            xnew=xi-alpha*gradx;
            ynew=yi-alpha*grady;
%          err = abs((ynew-yi)/yi);
      

            gvecx=[gvecx;xi];
            gvecy=[gvecy;yi];

            iter=iter+1;
       
            scatter(xi,yi,'g')
            if (abs(xnew-gx)<tol) && (abs(ynew-gy)<tol)
                disp('Goal Reached');
                disp(xi);
                disp(yi);
                break;
              
            end
            xi=xnew;
            yi=ynew;
            disp("x");
            disp(xi);
            disp("y");
            disp(yi);
          
%         disp("x is")
%         disp(x)
%         disp("y is")
%         disp(y)
          
%         end
%       end

      
end
hold off
 

%r=zeros(length(gvecx));
%plot3(gvecx,gvecy,r)
figure(3)
hold on 
%scatter(gvecx,gvecy);
quiver3(x,y,U,u,v,w,0.5)
contour(x,y,U);
circle(5,-3,sqrt(2)/2);
circle(6,6,sqrt(2)/2);
circle(gx,gx,.5)
circle(sx,sy,.5)
hold off

% plot3(gvecx,gvecy,U);
%hold off
% figure(2)
% hold on
% contour(x,y,U)
% quiver(x,y,Ux,Uy)
% hold off
%%
/
function circle (xc,yc,r)
ang=0:.01:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);
plot(xc+xp,yc+yp);
end
