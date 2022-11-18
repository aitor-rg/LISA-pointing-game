clear all;
close all;
rng(1,'twister');
MC = 1000;
dt = 0.1;
tspan = 0:dt:500;
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
v1 = randn(4,MC);
v2 = randn(4,MC);
v3 = randn(4,MC);
v4 = randn(4,MC);
len13 = 9*rand(4,MC);
len12 = (9-len13).*rand(4,MC);
y13v = v2./vecnorm(v2).*len13;
y12v = v1./vecnorm(v1).*len12;
y23v = (y12v-y13v)*0.8;
dphi1v = v4./vecnorm(v4).*0.8;
dphi1MC = 0.8*dphi1v;
dphi2MC = 0.8*dphi1v-y12v;
dphi3MC = 0.8*(dphi1v-y13v);

Umc = cell(MC); 
Xmc = cell(MC); 
Jmc = cell(MC); 
y12 = cell(MC);
y13 = cell(MC);
y23 = cell(MC);

parfor i = 1:MC
a = 0.5;
b = 0.2;

if mod(i/MC*100,1)==0
    disp([num2str(i/MC*100),'%']);
end
Rstep = 1;
Rsig = 1; 
gstep = 0.5*9;
gsig = 0.0645*9;
alpha =  0.93;
%% Sat position xyz (Don't matter in the code for now)
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



%% States and decision variables
c1 = s1;
c2 = s2;
c3 = s3;
rsphere = 1;


dphi1 = dphi1MC(:,i);
dphi2 = dphi2MC(:,i);
dphi3 = dphi3MC(:,i);

Wun = randn(3*7,1);
W = 1*sigma(1,Rsig,b,0.01*gsig)*Wun/norm(Wun);

U = zeros(3*7,length(tspan));
U(:,1) = [zeros(4,1);s1;zeros(4,1);s2;zeros(4,1);s3];

X = zeros(3*7,length(tspan));
X(:,1) = [dphi1;s1;dphi2;s2;dphi3;s3] + W;

M = zeros(3*7,length(tspan));


Uprev1 = zeros(4,1);
Uprev2 = zeros(4,1);
Uprev3 = zeros(4,1);
LL = 0;
J = zeros(3,length(tspan));
    x1 = X(1:7,1);
    x2 = X(8:14,1);
    x3 = X(15:21,1);
    J(1,1) = Cost1(x1,x2,x3,dphi1,dphi2,dphi3,(1/9^2));
    J(2,1) = Cost2(x1,x2,x3,dphi1,dphi2,dphi3,(1/9^2));
    J(3,1) = Cost3(x1,x2,x3,dphi1,dphi2,dphi3,(1/9^2));
bound1 = 0.3*(1-J(1,1));
bound2 = 0.3*(1-J(2,1));
bound3 = 0.3*(1-J(3,1));

for k = 2:length(tspan)

    x1 = X(1:7,k-1);
    x2 = X(8:14,k-1);
    x3 = X(15:21,k-1);

    J(1,k) = Cost1(x1,x2,x3,dphi1,dphi2,dphi3,(1/9^2));
    J(2,k) = Cost2(x1,x2,x3,dphi1,dphi2,dphi3,(1/9^2));
    J(3,k) = Cost3(x1,x2,x3,dphi1,dphi2,dphi3,(1/9^2));

    DJsmooth = 4*((J(1,k))-J(1,k-1))*W(1:4)/sigma(k-1-LL,Rsig,b,gsig);
    M(1:4,k) = DJsmooth;
    resf1 = (1-betaf)*resf1 + (betaf)*(M(1:4,k)-M(1:4,k-1));
    U(1:4,k) = U(1:4,k-1) + alpha*(U(1:4,k-1)-Uprev1) - stepsize(k-1-LL,Rstep,a,gstep)*(M(1:4,k)) ;

    DJsmooth = 4*(J(2,k)-J(2,k-1))*W(8:11)/sigma(k-1-LL,Rsig,b,gsig);
    M(8:11,k) = DJsmooth;
    resf2 = (1-betaf)*resf2 + (betaf)*(M(8:11,k)-M(8:11,k-1));
    U(8:11,k) = U(8:11,k-1) + alpha*(U(8:11,k-1)-Uprev2) - stepsize(k-1-LL,Rstep,a,gstep)*(M(8:11,k));

    DJsmooth = 4*(J(3,k)-J(3,k-1))*W(15:18)/sigma(k-1-LL,Rsig,b,gsig);
    M(15:18,k) = DJsmooth;
    U(15:18,k) = U(15:18,k-1) + alpha*(U(15:18,k-1)-Uprev3) - stepsize(k-1-LL,Rstep,a,gstep)*(M(15:18,k));
    Uprev1 = U(1:4,k-1);
    Uprev2 = U(8:11,k-1);
    Uprev3 = U(15:18,k-1);


    bound1 = bound1 + 0.01*(4.5-bound1);
    bound2 = bound2 + 0.01*(4.5-bound2);
    bound3 = bound3 + 0.01*(4.5-bound3);
    U(1:4,k) = min(max(U(1:4,k),-bound1),bound1);
    U(8:11,k) = min(max(U(8:11,k),-bound2),bound2);
    U(15:18,k) = min(max(U(15:18,k),-bound3),bound3);


    Wun = randn(3*7,1);
    W = Wun/norm(Wun);

    tau = -1/17;
    X(:,k) = exp(dt/tau)*X(:,k-1) + (1-exp(dt/tau))*(U(:,k)+sigma(k-LL,Rsig,b,gsig)*W);

end

Umc{i} = U; 
Xmc{i} = X; 
Jmc{i} = J; 
y12{i} = U(1:4,:)-U(8:11,:) -(-(dphi1MC(:,i)-dphi2MC(:,i))*ones(1,length(tspan)));
y13{i} = U(1:4,:)-U(15:18,:) -(-(dphi1MC(:,i)-dphi3MC(:,i))*ones(1,length(tspan)));
y23{i} = U(8:11,:)-U(15:18,:) -(-(dphi2MC(:,i)-dphi3MC(:,i))*ones(1,length(tspan)));
end
text = ['MC15nm',num2str(MC),'.mat'];
save(text,'Umc','Xmc','Jmc','dphi1MC','dphi2MC','dphi3MC','y12','y13','y23','tspan','-v7.3');
for i = 1:1000
y12m(i,:) = vecnorm(y12{i});
y13m(i,:) = vecnorm(y13{i});
y23m(i,:) = vecnorm(y23{i});
end
save('plotdataMC15t.mat','y12m','y13m','y23m','tspan','-v7.3');




