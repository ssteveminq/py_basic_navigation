% %Parameters
clc 
clear all 
close all
K=45; %attractive potential 
E=35; %repulsive potential
% area_width=12; %potential area width (m)

% K=50; %attractive potential 
% E=80; %repulsive potential
%>>>>>>> fc91092aa4cacc1af4cea39977739032d174c6de
sx=10; % start pos x
sy=11; %start pos y
gx=-1; %goal pos x
%<<<<<<< HEAD
gy=-5; %goal pos y
ox=[5.0 6.05]; %obstacle x pos
oy= [-3.0 6.05]; %obstacle y pos
% [x,y]= meshgrid(-12:.5:12); 
res=0.25;
[x,y] = meshgrid(-12:res:12,-12:res:12);
r=((x-gx).^2+(y-gy).^2).^.5;%Distance to Goal
%=======
%gy=-2; %goal pos y
%ox=[5 6]; %obstacle x pos
%oy= [-3 6]; %obstacle y pos
%[x,y]= meshgrid(-12:1:12); 
%[x,y] = meshgrid(-12:1.0:12,-12:1.0:12);
%r=((x-gx).^2+(y-gx).^2).^.5;%Distance to Goal
%>>>>>>> fc91092aa4cacc1af4cea39977739032d174c6de
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

iter=0;
% x=-2;
% y=2;

% maxiter=100;
errx=1;
erry=1;
%<<<<<<< HEAD
alpha=0.03;
tol=0.3;
%=======
%alpha=.25;
%tol=.75;
%>>>>>>> fc91092aa4cacc1af4cea39977739032d174c6de
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
% Udx=@(xs,ys)((50*(xs-1))/(xs^2-2*xs+ys^2+5-4*ys)^.5)-((80*(xs-5))/(xs^2-10*xs+ys^2+34-6*ys)^1.5)-((80*(xs-6))/(xs^2-12*xs+ys^2+72-12*ys)^1.5);
% Udy=@(xs,ys)((50*(ys-2))/(ys^2-4*ys+xs^2+5-2*xs)^.5)-((80*(ys-3))/(ys^2-6*ys+xs^2+34-10*xs)^1.5)-((80*(ys-6))/(ys^2-12*ys+xs^2+72-12*xs)^1.5);
% xs=10;
% ys=10;
% while xs~=gx && ys~=gy
%     xs=xs-alpha*Udx(xs,ys);
%     ys=ys-alpha*Udy(xs,ys);
%     iter=iter+1;
%     scatter(xs,ys,'g')
% end
% syms xs ys
% U=@(xs,ys)((K*((xs-gx)^2+(ys-gy)^2)^.5)+(E/((xs-xo1)^2+(ys-yo1)^2)^.5)+(E/((xs-xo2)^2+(ys-yo2)^2)^.5));

%<<<<<<< HEAD
xmin =-12.0;
ymin =-12.0;
maxiter=500;
%=======
%for i=1:25 %used to say -12:.5:12 
   
      %for j=1:25
%>>>>>>> fc91092aa4cacc1af4cea39977739032d174c6de

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
      

%             gvecx=[gvecx;xnew];
%             gvecy=[gvecy;ynew];

            iter=iter+1;
%<<<<<<< HEAD
       
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
%=======
%             %j=j+1;
%              while xi~=gx && yi~=gy
             %scatter(xi,yi,'g')
            %if xi<=gx+tol && xi>=gx-tol && yi<=gy+tol && yi>=gy-tol
                %disp('Goal Reached');
                %disp(xi);
                %disp(yi);
%                 
%                 
%         
            %end
            %xi=xnew;
            %yi=ynew;
%>>>>>>> fc91092aa4cacc1af4cea39977739032d174c6de
          
%         disp("x is")
%         disp(x)
%         disp("y is")
%         disp(y)
          
%<<<<<<< HEAD
%         end
%       end

      
%=======
%        end
%       end
%      
%       
%>>>>>>> fc91092aa4cacc1af4cea39977739032d174c6de
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
%<<<<<<< HEAD
%/
%=======
% syms xs ys
% Iter=0;
% Maxiter=25;
% Udx=((50*(xs-1))/(xs^2-2*xs+ys^2+5-4*ys)^.5)-((80*(xs-5))/(xs^2-10*xs+ys^2+34-6*ys)^1.5)-((80*(xs-6))/(xs^2-12*xs+ys^2+72-12*ys)^1.5);
% Udy=((50*(ys-2))/(ys^2-4*ys+xs^2+5-2*xs)^.5)-((80*(ys-3))/(ys^2-6*ys+xs^2+34-10*xs)^1.5)-((80*(ys-6))/(ys^2-12*ys+xs^2+72-12*xs)^1.5);
% while Iter<Maxiter
%     xs=xs-alpha*Udx(xs,yi);
%     ys=ys-alpha*Udx(xs,ys);
%     disp("Iter=");
%     disp(Iter);
%     disp("xi=")
%     disp(xs);
%     disp("yi=");
%     disp(ys);
% end
% %%
%>>>>>>> fc91092aa4cacc1af4cea39977739032d174c6de
function circle (xc,yc,r)
ang=0:.01:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);
plot(xc+xp,yc+yp);
end

