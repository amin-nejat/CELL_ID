function [S,R,T]=scaled_rotation(X,Y)
%Solves for Y = S*R*X + T
%where S is a diagonal scaling matrix
%R is a rotation matrix i.e. orthonormal and det(R)=1
%T is a translation

%% Remove NAN rows
idx=find(all(~isnan([X Y]),2));
X=X(idx,:);
Y=Y(idx,:);

%% MATLAB PROCRUSTES

% [~,~,transform]=procrustes(Y,X);
% S=transform.b;
% R=transform.T;
% T=transform.c(1,:);
% return

%% De-mean
Yhat=Y-mean(Y,1);
Xhat=X-mean(X,1);

%% Scale
sx = sqrt(sum(sum(Xhat.^2,2),1)/size(Xhat,1));
sy = sqrt(sum(sum(Yhat.^2,2),1)/size(Yhat,1));

% sx=sqrt(sum(Xhat.^2,1)/size(Xhat,1));
% sy=sqrt(sum(Yhat.^2,1)/size(Yhat,1));

Yhat=Yhat./sy;
Xhat=Xhat./sx;



%% Solve rotation
C = Yhat'*Xhat;
[U,~,V]=svd(C);

R0=V*U';

% Z=Xhat*R0;
% figure
% hold on
% plot3(Xhat(:,1),Xhat(:,2),Xhat(:,3),'o')
% plot3(Yhat(:,1),Yhat(:,2),Yhat(:,3),'o')
% plot3(Z(:,1),Z(:,2),Z(:,3),'o')
% grid on


%% Put it all together

S=diag(sy./sx);
R=R0;
T=(mean(Y,1)-mean(X,1)*R0*S);


%% Debug visual
% Z=X*R*S + T;
% figure(2)
% plot3(X(:,1),X(:,2),X(:,3),'o');grid on;
% hold on
% plot3(Y(:,1),Y(:,2),Y(:,3),'o');grid on;
% plot3(Z(:,1),Z(:,2),Z(:,3),'o');grid on;
end