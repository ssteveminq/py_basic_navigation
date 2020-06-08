[X,Y] = meshgrid(-12:.01:12);
Uatt= (.5.*((((X-6).^2+(Y+2).^2)).^2).*(((X-6).^2+(Y+2).^2).^.5)<=1)+((((X-6).^2+(Y+2).^2).^(1/2))-.5).*((((X-6).^2+(Y+2).^2).^.5)>1);
Urep=(5.*(((1./(((X).^2+(Y-1).^2).^.5))-(1/3)).*(((X).^2+(Y-1).^2).^.5)<=3)+(0).*((((X).^2+(Y-1).^2).^.5)>3));
U=Uatt-Urep;
%gamma=500;
%U= ((X-5).^2+(Y-6).^2)-(1./(((X-3).^2)+(((Y-4).^2)./4)-1));

% theta=2;% start x
% theta1=-2;%start  y
% alpha=.01;
% iter=0;
% errx=1;
% erry=1;
% tol=1e-6;
% maxiter=1000;
% Upx=diff(U,X);
% Upy=diff(U,Y);
% grad=[Upx,Upy];
% 
% [Upx,Upy]=gradient(U);
% figure;
% mesh(X,Y,Upx);



mesh(X,Y,U);
% hold on
% %contour(X,Y,Upx)
% hold off
% %quiver(X,Y,Upx,Upy);

% while iter<maxiter && errx>tol && erry>tol
%     temptheta=theta-alpha*diff(U,X);
%     temptheta1=theta-alpha*diff(U,Y);
%     errx=abs((theta-temptheta)/theta);
%     erry=abs((theta1-temptheta1)/theta1);
%     temptheta1=theta1;
%     temptheta=theta;
%     iter=iter+1;
% end
% if (iter<maxiter)
%     fprintf('X= %f and Y= %f after %d iterations \n',theta,theta1,iter)
% else 
%     fprintf('No solution after maxiter reached \n')
% end

