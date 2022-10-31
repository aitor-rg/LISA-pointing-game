clear all;
close all;
rng(0,'twister');
MC = 1000;
dt = 0.1;
tspan = 0:dt:5000;
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
v1 = randn(3,MC);
v2 = randn(3,MC);
v3 = randn(3,MC);
v4 = randn(3,MC);
len = rand(3,MC);
lenn = rand(3,MC);
y12v = v1./vecnorm(v1).*len*155;
y13v = v2./vecnorm(v2).*len*155;
y23v = v3./vecnorm(v3).*len*155;
dphi1v = v4./vecnorm(v4).*len*155;
Mat = [-eye(3) zeros(3);zeros(3) -eye(3); eye(3) -eye(3)];
Y = [y12v;y13v;y23v];
for i=1:MC
vec = Mat\(Y(:,i)-[dphi1v(:,i);dphi1v(:,i);zeros(3,1)]);
dphi2MC(:,i) = vec(1:3);
dphi3MC(:,i) = vec(4:6);
dphi1MC(:,i) = dphi1v(:,i);
end
% dphi1MC = 155*rand(3,MC);
% dphi2MC = 155*rand(3,MC);
% dphi3MC = 155*rand(3,MC);
Umc = cell(MC); 
Xmc = cell(MC); 
Jmc = cell(MC); 
y12 = cell(MC);
y13 = cell(MC);
y23 = cell(MC);

for i = 1:MC
%% gradient parameters a = 0.5, b = 0.2, R = 1E4
a = 0.5; %a=0.6
b = 0.3; % b = 0.07; a+2b in (0.5,1]
% a+2*b
if mod(MC,i)==0
    disp([num2str(i/MC*100),'%']);
end
Rstep = 1;
Rsig = 1; % 1E-7 also works
gstep = 200;
gsig = 20;
% smoothing parameter
beta = 0;
alpha = 0.95;
%% Sat position xyz
L = 100;
l = [L,0,0]';
s1 = Rz(90*pi/180)*l;
s2 = Rz(-30*pi/180)*l;
s3 = Rz(210*pi/180)*l;

%% Euler angles references
% absolutes
r_phi1 = [0,0,-90*pi/180]';
r_phi2 = [0,0,150*pi/180]';
r_phi3 = [0,0,30*pi/180]';

% relatives
r_phi12 = [0,0,(150-(-90))*pi/180];
r_phi13 = [0,0,(30-(-90))*pi/180];

r_phi21 = [0,0,(-90-150)*pi/180];
r_phi23 = [0,0,(30-150)*pi/180];

r_phi31 = [0,0,(-90-30)*pi/180];
r_phi32 = [0,0,(150-30)*pi/180];

% Rot matrices (ref frames)
% S1 = [Rz(phi1(3)),s1; zeros(1,3), 1];
% S2 = [Rz(phi2(3)),s2; zeros(1,3), 1];
% S3 = [Rz(phi3(3)),s3; zeros(1,3), 1];

% hold on;
% triangle = [s1';s2';s3';s1'];
% plot3(triangle(:,1),triangle(:,2),triangle(:,3));
% trplot(S1,'length',30);
% trplot(S2,'length',30);
% trplot(S3,'length',30);

%% States and decision variables
c1 = s1;
c2 = s2;
c3 = s3;
rsphere = 1;

% Positions
s1 = s1 + 0*randn(3,1);
s2 = s2 + 0*randn(3,1);
s3 = s3 + 0*randn(3,1);
% Relative Euler angles IN MICRORADIANTS?
dphi1 = dphi1MC(:,i);
dphi2 = dphi2MC(:,i);
dphi3 = dphi3MC(:,i);
% Noise
Wun = randn(3*6,1);
W = sigma(1,Rsig,b,gsig)*Wun/norm(Wun);

U = zeros(3*6,length(tspan));
U(:,1) = [zeros(3,1);s1;zeros(3,1);s2;zeros(3,1);s3];

X = zeros(3*6,length(tspan));
X(:,1) = [dphi1;s1;dphi2;s2;dphi3;s3] + W;

M = zeros(3*6,length(tspan));

% U = zeros(6*3,length(tspan));
% U(:,1) = [phi1; s1; phi2; s2; phi3; s3];
% X = zeros(6*3,length(tspan));
% X(:,1) = [phi1; s1; phi2; s2; phi3; s3] + W;
Uprev1 = zeros(3,1);
Uprev2 = zeros(3,1);
Uprev3 = zeros(3,1);
LL = 0;
J = zeros(3,length(tspan));
res = 0;
gamma = 1; 
for k = 2:length(tspan)
%     x1 = X(1:6,k-1);   %s1
%     x2 = X(7:12,k-1);  %s2
%     x3 = X(13:18,k-1); %s3
    x1 = X(1:6,k-1);
    x2 = X(7:12,k-1);
    x3 = X(13:18,k-1);

    J(1,k) = Cost1(x1,x2,x3,dphi1,dphi2,dphi3,(1/155^2));
    J(2,k) = Cost2(x1,x2,x3,dphi1,dphi2,dphi3,(1/155^2));
    J(3,k) = Cost3(x1,x2,x3,dphi1,dphi2,dphi3,(1/155^2));

    
    % Min
    DJsmooth = 3*gamma*((J(1,k)))*W(1:3)/sigma(k-1-LL,Rsig,b,gsig);
    M(1:3,k) = beta*M(1:3,k-1) + DJsmooth*(1-beta);
    U(1:3,k) = (U(1:3,k-1) - stepsize(k-1-LL,Rstep,a,gstep)*(M(1:3,k)) + (0/k^0)*(U(1:3,k-1)-Uprev1))/1;

    DJsmooth = 3*gamma*(J(2,k))*W(7:9)/sigma(k-1-LL,Rsig,b,gsig);
    M(7:9,k) = beta*M(7:9,k-1) + DJsmooth*(1-beta);
    U(7:9,k) = (U(7:9,k-1) - stepsize(k-1-LL,Rstep,a,gstep)*(M(7:9,k)) + (0/k^0)*(U(7:9,k-1)-Uprev2))/1;

    DJsmooth = 3*gamma*(J(3,k))*W(13:15)/sigma(k-1-LL,Rsig,b,gsig);
    M(13:15,k) = beta*M(13:15,k-1) + DJsmooth*(1-beta);
    U(13:15,k) = (U(13:15,k-1) - stepsize(k-1-LL,Rstep,a,gstep)*(M(13:15,k)) + (0/k^0)*(U(13:15,k-1)-Uprev3))/1;
    Uprev1 = U(1:3,k-1);
    Uprev2 = U(7:9,k-1);
    Uprev3 = U(13:15,k-1);
    U(:,k) = min(max(U(:,k),-80),80);
%     voffset = 10*1e9;
%     U(1:3,k)   = max(min(U(1:3,k),phi1*1e9+voffset),phi1*1e9-voffset);
%     U(7:9,k)   = max(min(U(7:9,k),phi2*1e9+voffset),phi2*1e9-voffset);
%     U(13:15,k) = max(min(U(13:15,k),phi3*1e9+voffset),phi3*1e9-voffset);
%
%     U(4:6,k) = c1 + rsphere*(U(4:6,k)-c1)/norm((U(4:6,k)-c1));
%     U(10:12,k) = c2 + rsphere*(U(10:12,k)-c2)/norm((U(10:12,k)-c2));
%     U(16:18,k) = c3 + rsphere*(U(16:18,k)-c3)/norm((U(16:18,k)-c3));

    Wun = randn(3*6,1);
    W = Wun/norm(Wun);
%     [~,q] = ode45(@(t,x) model(t,x,U(:,k)+W),[tspan(k-1),tspan(k)],X(:,k-1),opts);
%     X(:,k) = q(end,:)';
    tau = -0.01;
    X(:,k) = exp(dt/tau)*X(:,k-1) + (1-exp(dt/tau))*(U(:,k)+sigma(k-LL,Rsig,b,gsig)*W);
%     if(norm([U(1,k)-U(7,k)+(dphi1(1)-dphi2(1));U(2,k)-U(8,k)+(dphi1(2)-dphi2(2));U(3,k)-U(9,k)+(dphi1(3)-dphi2(3))])<10^(-(9+res)))
%         LL = 0;round(k/3);
%         res = 1;
%         %gamma = 100*gamma;
% %         Uprev1 = gamma*U(1:3,k-1);
% %         Uprev2 = gamma*U(7:9,k-1);
% %         Uprev3 = gamma*U(13:15,k-1);
%     end
end

% figure(1);
% subplot(3,1,1);
% plot(J(1,2:end));
% subplot(3,1,2);
% plot(J(2,2:end));
% subplot(3,1,3);
% plot(J(3,2:end));
% grid on;

% figure(1);
% subplot(3,1,1);
% plot(M(1,2:end));
% subplot(3,1,2);
% plot(M(2,2:end));
% subplot(3,1,3);
% plot(M(3,2:end));
% grid on;
% 
% figure(2);
% subplot(3,1,1);
% title('Errors')
% plot(U(7,:)-U(1,:));
% subplot(3,1,2);
% plot(U(8,:)-U(2,:));
% subplot(3,1,3);
% plot((150-(-90))*pi/180*1e9 -(U(9,:)-U(3,:)));
% grid on;
% figure(3);
% subplot(3,1,1);
% plot(U(13,:)-U(1,:));
% subplot(3,1,2);
% plot(U(14,:)-U(2,:));
% subplot(3,1,3);
% plot((30-(-90))*pi/180*1e9 -(U(15,:)-U(3,:)));
% grid on;
% figure(4);
% subplot(3,1,1);
% plot(U(13,:)-U(7,:));
% subplot(3,1,2);
% plot(U(14,:)-U(8,:));
% subplot(3,1,3);
% plot((30-150)*pi/180*1e9 -(U(15,:)-U(9,:)));
% grid on;


% figure(3);
% subplot(3,1,1);
% title('Angles')
% plot(X(1,:)-X(7,:));hold on;
% plot(U(1,:)-U(7,:));hold on;
% plot(-(dphi1(1)-dphi2(1))*ones(1,length(tspan)));
% subplot(3,1,2);
% plot(X(2,:)-X(8,:));hold on;
% plot(U(2,:)-U(8,:));hold on;
% plot(-(dphi1(2)-dphi2(2))*ones(1,length(tspan)));
% subplot(3,1,3);
% plot(X(3,:)-X(9,:));hold on;
% plot(U(3,:)-U(9,:));hold on;
% plot(-(dphi1(3)-dphi2(3))*ones(1,length(tspan)));
% grid on;
% 
% figure(4);
% subplot(3,1,1);
% title('Angles')
% plot(X(1,:)-X(13,:));hold on;
% plot(U(1,:)-U(13,:));hold on;
% plot(-(dphi1(1)-dphi3(1))*ones(1,length(tspan)));
% subplot(3,1,2);
% plot(X(2,:)-X(14,:));hold on;
% plot(U(2,:)-U(14,:));hold on;
% plot(-(dphi1(2)-dphi3(2))*ones(1,length(tspan)));
% subplot(3,1,3);
% plot(X(3,:)-X(15,:));hold on;
% plot(U(3,:)-U(15,:));hold on;
% plot(-(dphi1(3)-dphi3(3))*ones(1,length(tspan)));
% grid on;
% 
% figure(5);
% subplot(3,1,1);
% title('Angles')
% plot(X(7,:)-X(13,:));hold on;
% plot(U(7,:)-U(13,:));hold on;
% plot(-(dphi2(1)-dphi3(1))*ones(1,length(tspan)));
% subplot(3,1,2);
% plot(X(8,:)-X(14,:));hold on;
% plot(U(8,:)-U(14,:));hold on;
% plot(-(dphi2(2)-dphi3(2))*ones(1,length(tspan)));
% subplot(3,1,3);
% plot(X(9,:)-X(15,:));hold on;
% plot(U(9,:)-U(15,:));hold on;
% plot(-(dphi2(3)-dphi3(3))*ones(1,length(tspan)));
% grid on;
% 
% figure(6);
% subplot(3,1,1);
% title('Angles')
% plot(U(1,:)-U(7,:)-(-(dphi1(1)-dphi2(1))*ones(1,length(tspan))));hold on;
% grid on
% subplot(3,1,2);
% plot(U(2,:)-U(8,:)-(-(dphi1(2)-dphi2(2))*ones(1,length(tspan))));hold on;
% grid on 
% subplot(3,1,3);
% plot(U(3,:)-U(9,:)-(-(dphi1(3)-dphi2(3))*ones(1,length(tspan))));hold on;
% grid on;
Umc{i} = U; 
Xmc{i} = X; 
Jmc{i} = J; 
y12{i} = X(1:3,:)-X(7:9,:) -(-(dphi1MC(:,i)-dphi2MC(:,i))*ones(1,length(tspan)));
y13{i} = X(1:3,:)-X(13:15,:) -(-(dphi1MC(:,i)-dphi3MC(:,i))*ones(1,length(tspan)));
y23{i} = X(7:9,:)-X(13:15,:) -(-(dphi2MC(:,i)-dphi3MC(:,i))*ones(1,length(tspan)));
end
text = ['MC',num2str(MC),'nomomentum','.mat'];
save(text,'Umc','Xmc','Jmc','dphi1MC','dphi2MC','dphi3MC','y12','y13','y23','tspan','-v7.3');
%% differential equation
function dy = model(t,x,u)
    tau = 0.01;
    y = u-x;
    dy = (1/tau)*y;
end




