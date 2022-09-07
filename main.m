clear all;
close all;
rng(1);

dt = 0.001;
tspan = 0:dt:300;
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);

%% gradient parameters a = 0.5, b = 0.2, R = 1E4
a = 0.6;
b = 0.07; % a+2b in (0.5,1]
a+2*b
Rstep = 1E-2;
Rsig = 1E-2; % 1E-7 also works
% smoothing parameter
beta = 0.7;

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
dphi1 = 1*randn(3,1);
dphi2 = 1*randn(3,1);
dphi3 = 1*randn(3,1);
% Noise
W = sigma(1,Rsig,b)*randn(3*6,1);

U = zeros(3*6,length(tspan));
U(:,1) = [zeros(3,1);s1;zeros(3,1);s2;zeros(3,1);s3];

X = zeros(3*6,length(tspan));
X(:,1) = [dphi1;s1;dphi2;s2;dphi3;s3] + W;

M = zeros(3*6,length(tspan));

% U = zeros(6*3,length(tspan));
% U(:,1) = [phi1; s1; phi2; s2; phi3; s3];
% X = zeros(6*3,length(tspan));
% X(:,1) = [phi1; s1; phi2; s2; phi3; s3] + W;

J = zeros(3,length(tspan));
for k = 2:length(tspan)

%     x1 = X(1:6,k-1);   %s1
%     x2 = X(7:12,k-1);  %s2
%     x3 = X(13:18,k-1); %s3
    x1 = X(1:6,k-1);
    x2 = X(7:12,k-1);
    x3 = X(13:18,k-1);

    J(1,k) = Cost1(x1,x2,x3,dphi1,dphi2,dphi3);
    J(2,k) = Cost2(x1,x2,x3,dphi1,dphi2,dphi3);
    J(3,k) = Cost3(x1,x2,x3,dphi1,dphi2,dphi3);

    
    % Min
    DJsmooth = J(1,k)*W(1:3)*sigma(k,Rsig,b)^2/sigma(k-1,Rsig,b)^2;
    M(1:3,k) = beta*M(1:3,k-1) + DJsmooth*(1-beta);
    U(1:3,k) = U(1:3,k-1) + stepsize(k,Rstep,a)*M(1:3,k);

    DJsmooth = J(2,k)*W(7:9)*sigma(k,Rsig,b)^2/sigma(k-1,Rsig,b)^2;
    M(7:9,k) = beta*M(7:9,k-1) + DJsmooth*(1-beta);
    U(7:9,k) = U(7:9,k-1) + stepsize(k,Rstep,a)*M(7:9,k);

    DJsmooth = J(3,k)*W(13:15)*sigma(k,Rsig,b)^2/sigma(k-1,Rsig,b)^2;
    M(13:15,k) = beta*M(13:15,k-1) + DJsmooth*(1-beta);
    U(13:15,k) = U(13:15,k-1) + stepsize(k,Rstep,a)*M(13:15,k);
    
%     voffset = 10*1e9;
%     U(1:3,k)   = max(min(U(1:3,k),phi1*1e9+voffset),phi1*1e9-voffset);
%     U(7:9,k)   = max(min(U(7:9,k),phi2*1e9+voffset),phi2*1e9-voffset);
%     U(13:15,k) = max(min(U(13:15,k),phi3*1e9+voffset),phi3*1e9-voffset);
%
%     U(4:6,k) = c1 + rsphere*(U(4:6,k)-c1)/norm((U(4:6,k)-c1));
%     U(10:12,k) = c2 + rsphere*(U(10:12,k)-c2)/norm((U(10:12,k)-c2));
%     U(16:18,k) = c3 + rsphere*(U(16:18,k)-c3)/norm((U(16:18,k)-c3));

    W = sigma(k,Rsig,b)*randn(3*6,1);
    %[~,q] = ode45(@(t,x) model(t,x,U(:,k)+W),[tspan(k-1),tspan(k)],X(:,k-1),opts);
    X(:,k) = U(:,k)+W;%q(end,:)';
end

figure(1);
subplot(3,1,1);
plot(J(1,2:end));
subplot(3,1,2);
plot(J(2,2:end));
subplot(3,1,3);
plot(J(3,2:end));
grid on;

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


figure(3);
subplot(3,1,1);
title('Angles')
plot(X(1,:)-X(7,:));hold on;
plot(U(1,:)-U(7,:));hold on;
plot(-(dphi1(1)-dphi2(1))*ones(1,length(tspan)));
subplot(3,1,2);
plot(X(2,:)-X(8,:));hold on;
plot(U(2,:)-U(8,:));hold on;
plot(-(dphi1(2)-dphi2(2))*ones(1,length(tspan)));
subplot(3,1,3);
plot(X(3,:)-X(9,:));hold on;
plot(U(3,:)-U(9,:));hold on;
plot(-(dphi1(3)-dphi2(3))*ones(1,length(tspan)));
grid on;
% subplot(3,1,2);
% plot(U(7,:));hold on;
% plot(U(8,:));
% plot(U(9,:));
% grid on;
% subplot(3,1,3);
% plot(U(13,:));hold on;
% plot(U(14,:));
% plot(U(15,:));
% grid on;


% figure(4);
% subplot(3,1,1);
% title('Positions')
% plot(U(4,:));hold on;
% plot(U(5,:));
% plot(U(6,:));
% grid on;
% subplot(3,1,2);
% plot(U(10,:));hold on;
% plot(U(11,:));
% plot(U(12,:));
% grid on;
% subplot(3,1,3);
% plot(U(16,:));hold on;
% plot(U(17,:));
% plot(U(18,:));
% grid on;



% figure(1)
% triangle = [X(4:6,end)';X(10:12,end)';X(16:18,end)';X(4:6,end)'];
% plot3(triangle(:,1),triangle(:,2),triangle(:,3),'r');
% R1 = Rz(X(3,end)*1e-9)*Ry(X(2,end)*1e-9)*Rx(X(1,end)*1e-9);
% R2 = Rz(X(9,end)*1e-9)*Ry(X(8,end)*1e-9)*Rx(X(7,end)*1e-9);
% R3 = Rz(X(15,end)*1e-9)*Ry(X(14,end)*1e-9)*Rx(X(13,end)*1e-9);
% T1 = [R1,X(4:6,end); zeros(1,3), 1];
% T2 = [R2,X(10:12,end); zeros(1,3), 1];
% T3 = [R3,X(16:18,end); zeros(1,3), 1];
% trplot(T1,'length',30,'color','r');
% trplot(T2,'length',30,'color','r');
% trplot(T3,'length',30,'color','r');
% axis equal


%figure
% for i =1:length(U(4,:))
% scatter3(U(4,i),U(5,i),U(6,i),100, 'LineWidth', 3); hold on;
% scatter3(U(10,i),U(11,i),U(12,i),100, 'LineWidth', 3); hold on;
% scatter3(U(16,i),U(17,i),U(18,i),100, 'LineWidth', 3); hold on;
% axis equal
% pause(0.01)
% end

%% differential equation
function dy = model(t,x,u)
    tau = 0.1;
    y = u-x;
    dy = (1/tau)*y;
end




