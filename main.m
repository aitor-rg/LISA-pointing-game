clear all;
close all;
rng(0);

l = 100;
dt = 0.1;
tspan = 0:dt:500;
opts = odeset('RelTol',1e-7,'AbsTol',1e-7);

% gradient parameters
a = 0.55;
b = 0.2; % a+2b in (0.5,1]
a+2*b
Rstep = 1E4;
Rsig = 1E4;
% smoothing parameter
beta = 0.99; 

% Sat position xyz
s1 = [l,0,0]' + 0.1*randn(3,1);
s2 = Rz(60*pi/180)*s1 + 0.1*randn(3,1);
s3 = Rz(60*pi/180)*s2 + 0.1*randn(3,1);
% Euler angles
phi1 = [0,0,-90*pi/180]' + 0.1*randn(3,1);
phi2 = [0,0,120*pi/180]' + 0.1*randn(3,1);
phi3 = [0,0,30*pi/180]' + 0.1*randn(3,1);
% Rot matrices (ref frames)
S1 = Rz(phi1(3));
S2 = Rz(phi2(3));
S3 = Rz(phi3(3));

W = sigma(1,Rsig,b)*randn(9,1);

U = zeros(9,length(tspan));
U(:,1) = [phi1;phi2;phi3];

X = zeros(9,length(tspan));
X(:,1) = U(:,1) + W;

M = zeros(9,length(tspan));

% U = zeros(6*3,length(tspan));
% U(:,1) = [phi1; s1; phi2; s2; phi3; s3];
% X = zeros(6*3,length(tspan));
% X(:,1) = [phi1; s1; phi2; s2; phi3; s3] + W;

J = zeros(3,length(tspan));
for k = 2:length(tspan)    
    
%     x1 = X(1:6,k-1);   %s1
%     x2 = X(7:12,k-1);  %s2
%     x3 = X(13:18,k-1); %s3
    x1 = X(1:3,k-1);
    x2 = X(4:6,k-1);
    x3 = X(7:9,k-1);
    
    J(1,k) = Cost1(x1,x2,x3);
    J(2,k) = Cost2(x1,x2,x3);
    J(3,k) = Cost3(x1,x2,x3);
    
    % Min
    DJsmooth = J(1,k)*W(1:3)*sigma(k,Rsig,b)^2/sigma(k-1,Rsig,b)^2;
    M(1:3,k) = beta*M(1:3,k-1) + DJsmooth*(1-beta);
    U(1:3,k) = U(1:3,k-1) - stepsize(k,Rstep,a)*M(1:3,k);
    
    DJsmooth = J(2,k)*W(4:6)*sigma(k,Rsig,b)^2/sigma(k-1,Rsig,b)^2;
    M(4:6,k) = beta*M(4:6,k-1) + DJsmooth*(1-beta);
    U(4:6,k) = U(4:6,k-1) - stepsize(k,Rstep,a)*M(4:6,k);
    
    DJsmooth = J(3,k)*W(7:9)*sigma(k,Rsig,b)^2/sigma(k-1,Rsig,b)^2;
    M(7:9,k) = beta*M(7:9,k-1) + DJsmooth*(1-beta);
    U(7:9,k) = U(7:9,k-1) - stepsize(k,Rstep,a)*M(7:9,k);
    
    U(1:3,k)   = max(min(U(1:3,k),pi),-pi);
    U(4:6,k)   = max(min(U(4:6,k),pi),-pi);
    U(7:9,k) = max(min(U(7:9,k),pi),-pi);
    
    
    W = sigma(k,Rsig,b)*randn(9,1);
    [~,q] = ode45(@(t,x) model(t,x,U(:,k)+W),[tspan(k-1),tspan(k)],X(:,k-1),opts);
    X(:,k) = q(end,:)';
end

figure;
subplot(3,1,1);
plot(J(1,:));
subplot(3,1,2);
plot(J(2,:));
subplot(3,1,3);
plot(J(3,:));
grid on;

figure;
subplot(3,1,1);
plot(U(1,:)*180/pi);hold on;
plot(U(2,:)*180/pi);
plot(U(3,:)*180/pi);
grid on;
subplot(3,1,2);
plot(U(4,:)*180/pi);hold on;
plot(U(5,:)*180/pi);
plot(U(6,:)*180/pi);
grid on;
subplot(3,1,3);
plot(U(7,:)*180/pi);hold on;
plot(U(8,:)*180/pi);
plot(U(9,:)*180/pi);
grid on;


%% differential equation
function dy = model(t,x,u)
    tau = 0.0001;
    y = u-x;
    dy = (1/tau)*y;
end




