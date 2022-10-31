function J = Const23(s1,s2,s3,dphi1,dphi2,dphi3,sig)
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
    
    %nr = Rz(30*pi/180)*Rz(phi1(3)*1e-9)*Ry(phi1(2)*1e-9)*Rx(phi1(1)*1e-9)*[1;0;0];
    
    %nl = Rz(-30*pi/180)*Rz(phi1(3)*1e-9)*Ry(phi1(2)*1e-9)*Rx(phi1(1)*1e-9)*[1;0;0];
    
    %laser direction of s2
    %dr = Rz(-30*pi/180)*Rz(phi2(3)*1e-9)*Ry(phi2(2)*1e-9)*Rx(phi2(1)*1e-9)*[1;0;0];
    %laser direction of s3
    %dl = Rz(30*pi/180)*Rz(phi3(3)*1e-9)*Ry(phi3(2)*1e-9)*Rx(phi3(1)*1e-9)*[1;0;0];
    %offset_r = xyz2-xyz1;
    %offset_l = xyz3-xyz1;
    %tr = dot(nr,offset_r)/dot(nr,dr);
    %tl = dot(nl,offset_l)/dot(nl,dl);
    

    siga = sig; %inv(50*eye(3));
    %sigp = 10;

    %Jr = 1*exp(-0.5*(offset_r - dr*tr)'*sig^-1*(offset_r - dr*tr))-1*(dot(nr,dr)+1)^2;
    %Jl = 1*exp(-0.5*(offset_l - dl*tl)'*sig^-1*(offset_l - dl*tl))-1*(dot(nl,dl)+1)^2;

    x23 = phi2-phi3;
    r23 = -(dphi2-dphi3);
%     if (x12-r12)<1e-2
%         factor = 1;
%     else
%         factor = 1;
%     end
    Jr = ((x23-r23))'*((x23-r23));

    J = Jr - sig^2;
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