function J = Cost1(s1,s2,s3,dphi1,dphi2,dphi3)
    % s1
    phi1 = s1(1:3);
    xyz1 = s1(4:6);
    % left of s1
    phi2 = s2(1:3);
    xyz2 = s2(4:6);
    % right of s1
    phi3 = s3(1:3);
    xyz3 = s3(4:6);
   
    % panels of s1
    
    nr = Rz(30*pi/180)*Rz(phi1(3)*1e-9)*Ry(phi1(2)*1e-9)*Rx(phi1(1)*1e-9)*[1;0;0];
    
    nl = Rz(-30*pi/180)*Rz(phi1(3)*1e-9)*Ry(phi1(2)*1e-9)*Rx(phi1(1)*1e-9)*[1;0;0];
    
    %laser direction of s2
    dr = Rz(-30*pi/180)*Rz(phi2(3)*1e-9)*Ry(phi2(2)*1e-9)*Rx(phi2(1)*1e-9)*[1;0;0];
    %laser direction of s3
    dl = Rz(30*pi/180)*Rz(phi3(3)*1e-9)*Ry(phi3(2)*1e-9)*Rx(phi3(1)*1e-9)*[1;0;0];
    offset_r = xyz2-xyz1;
    offset_l = xyz3-xyz1;
    tr = dot(nr,offset_r)/dot(nr,dr);
    tl = dot(nl,offset_l)/dot(nl,dl);
    

    siga = inv(50*eye(3));
    sigp = 10;

    %Jr = 1*exp(-0.5*(offset_r - dr*tr)'*sig^-1*(offset_r - dr*tr))-1*(dot(nr,dr)+1)^2;
    %Jl = 1*exp(-0.5*(offset_l - dl*tl)'*sig^-1*(offset_l - dl*tl))-1*(dot(nl,dl)+1)^2;

    x12 = phi1-phi2;
    r12 = -(dphi1-dphi2);
    x13 = phi1-phi3;
    r13 = -(dphi1-dphi3);
     
    Jr = -0.5*(x12-r12)'*siga*(x12-r12);
    Jl = -0.5*(x13-r13)'*siga*(x13-r13);

    J = Jr + Jl;
%     r1 = [0,0,-90*pi/180]';
%     
%     c12 = [0,0,210*pi/180]';
%     c13 = [0,0,120*pi/180]';
%     
%     J_12 = (phi2-phi1-c12)'*(phi2-phi1-c12);
%     J_13 = (phi3-phi1-c13)'*(phi3-phi1-c13);
%     
%     J = (J_12 + J_13) + 1*(phi1-r1)'*(phi1-r1);
end