function J = Cost2(s1,s2,s3,dphi1,dphi2,dphi3)
    % right of s2
    phi1 = s1(1:3);
    xyz1 = s1(4:6);
    % s2
    phi2 = s2(1:3);
    xyz2 = s2(4:6);
    % left of s2
    phi3 = s3(1:3);
    xyz3 = s3(4:6);
    
    % panels of s2
    nr = Rz(-30*pi/180)*Rz(phi2(3)*1e-9)*Ry(phi2(2)*1e-9)*Rx(phi2(1)*1e-9)*[1;0;0];
    nl = Rz(30*pi/180)*Rz(phi2(3)*1e-9)*Ry(phi2(2)*1e-9)*Rx(phi2(1)*1e-9)*[1;0;0];
    
    % laser direction of s1
    dr = Rz(30*pi/180)*Rz(phi1(3)*1e-9)*Ry(phi1(2)*1e-9)*Rx(phi1(1)*1e-9)*[1;0;0];
    % laser direction of s3
    dl = Rz(-30*pi/180)*Rz(phi3(3)*1e-9)*Ry(phi3(2)*1e-9)*Rx(phi3(1)*1e-9)*[1;0;0];
    
    offset_r = xyz1-xyz2;
    offset_l = xyz3-xyz2;
    
    tr = dot(nr,offset_r)/dot(nr,dr);
    tl = dot(nl,offset_l)/dot(nl,dl);
    
    %Jr = norm(offset_r - dr*tr)^2 - dot(nr,dr)^2;
    %Jl = norm(offset_l - dl*tl)^2 - dot(nl,dl)^2;
    siga = inv(50*eye(3));
    sigp = 10;
    
    x21 = phi2-phi1;
    r21 = -(dphi2-dphi1);
    x23 = phi2-phi3;
    r23 = -(dphi2-dphi3);
     
    Jr = -0.5*(x21-r21)'*siga*(x21-r21);
    Jl = -0.5*(x23-r23)'*siga*(x23-r23);
    
    J = Jr + Jl ;

%     r2 = [0,0,120*pi/180]';
%     
%     c21 = [0,0,-210*pi/180]';
%     c23 = [0,0,-90*pi/180]';
%     
%     J_21 = (phi1-phi2-c21)'*(phi1-phi2-c21);
%     J_23 = (phi3-phi2-c23)'*(phi3-phi2-c23);
%     
%     J = (J_21 + J_23) + 1*(phi2-r2)'*(phi2-r2);

end