l = 100;
dt = 0.1;
sigma = 0.0001;
stepsize = 0.0001;
tspan = 0:dt:100;
opts = odeset('RelTol',1e-7,'AbsTol',1e-7);

% Sat position xyz
s1 = [l,0,0]' + 0.1*randn(3,1);
s2 = Rz(60*pi/180)*s1 + 0.1*randn(3,1);
s3 = Rz(60*pi/180)*s2 + 0.1*randn(3,1);
% Euler angles
phi1 = [0,0,-pi/2]' + 0.001*randn(3,1);
phi2 = [0,0,120*pi/180]' + 0.001*randn(3,1);
phi3 = [0,0,-120*pi/180]' + 0.001*randn(3,1);
% Rot matrices (ref frames)
S1 = Rz(phi1(3));
S2 = Rz(phi2(3));
S3 = Rz(phi3(3));

W = sigma*randn(18,1);

U = zeros(6*3,length(tspan));
U(:,1) = [phi1; s1; phi2; s2; phi3; s3];
X = zeros(6*3,length(tspan));
X(:,1) = [phi1; s1; phi2; s2; phi3; s3] + W;

J = zeros(3,length(tspan));
for k = 2:length(tspan)    
    
    x1 = X(1:6,k-1);   %s1
    x2 = X(7:12,k-1);  %s2
    x3 = X(13:18,k-1); %s3
    
    J(1,k) = Cost1(x1,x2,x3);
    J(2,k) = Cost2(x1,x2,x3);
    J(3,k) = Cost3(x1,x2,x3);
    
    % Min
    U(1:6,k)   = U(1:6,k-1)   - J(1,k)*W(1:6)*stepsize;
    U(7:12,k)  = U(7:12,k-1)  - J(2,k)*W(7:12)*stepsize;
    U(13:18,k) = U(13:18,k-1) - J(3,k)*W(13:18)*stepsize;

    U(1:3,k)   = max(min(U(1:3,k),pi),-pi);
    U(7:9,k)   = max(min(U(7:9,k),pi),-pi);
    U(13:15,k) = max(min(U(13:15,k),pi),-pi);
    
    W = sigma*randn(18,1);
    [~,q] = ode45(@(t,x) model(t,x,U(:,k)+W),[tspan(k-1),tspan(k)],X(:,k-1),opts);
    X(:,k) = q(end,:)';
end

plot(J(1,:));hold on;
plot(J(2,:));
plot(J(3,:));

%% differential equation
function dy = model(t,x,u)
    tau = 0.001;
    y = u-x;
    dy = (1/tau)*y;
end




