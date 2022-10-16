function J = Cost3(s1,s2,s3,dphi1,dphi2,dphi3,sig)
    % right of s2
    phi1 = s1(1:3);
    xyz1 = s1(4:6);
    % s2
    phi2 = s2(1:3);
    xyz2 = s2(4:6);
    % s3
    phi3 = s3(1:3);
    xyz3 = s3(4:6);
    
    % panels of s2
    nr = Rz(-30*pi/180)*Rz(phi3(3)*1e-9)*Ry(phi3(2)*1e-9)*Rx(phi3(1)*1e-9)*[1;0;0];
    nl = Rz(30*pi/180)*Rz(phi3(3)*1e-9)*Ry(phi3(2)*1e-9)*Rx(phi3(1)*1e-9)*[1;0;0];
    
    % laser direction of s1
    dr = Rz(30*pi/180)*Rz(phi2(3)*1e-9)*Ry(phi2(2)*1e-9)*Rx(phi2(1)*1e-9)*[1;0;0];
    % laser direction of s3
    dl = Rz(-30*pi/180)*Rz(phi1(3)*1e-9)*Ry(phi1(2)*1e-9)*Rx(phi1(1)*1e-9)*[1;0;0];
    
    offset_r = xyz2-xyz3;
    offset_l = xyz1-xyz3;
    
    tr = dot(nr,offset_r)/dot(nr,dr);
    tl = dot(nl,offset_l)/dot(nl,dl);
    
    %Jr = norm(offset_r - dr*tr)^2 - dot(nr,dr)^2;
    %Jl = norm(offset_l - dl*tl)^2 - dot(nl,dl)^2;
   
    siga = sig;%inv(50*eye(3));
    sigp = 10;
    
    x31 = phi3-phi1;
    r31 = -(dphi3-dphi1);
    x32 = phi3-phi2;
    r32 = -(dphi3-dphi2);
     
    Jr = 0.5*(x31-r31)'*siga*(x31-r31);
    Jl = 0.5*(x32-r32)'*siga*(x32-r32);
    
    J = Jr + Jl ;

%     r3 = [0,0,30*pi/180]';
%  
%     c31 = [0,0,-120*pi/180]';
%     c32 = [0,0,90*pi/180]';
%     
%     J_31 = (phi1-phi3-c31)'*(phi1-phi3-c31);
%     J_32 = (phi2-phi3-c32)'*(phi2-phi3-c32);
%     
%     J = (J_31 + J_32) + 1*(phi3-r3)'*(phi3-r3);
end