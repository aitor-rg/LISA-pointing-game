function J = Cost2(s1,s2,s3,dphi1,dphi2,dphi3,sig)
    % right of s2
    phi1 = s1(1:4);
    xyz1 = s1(5:7);
    % left of s1
    phi2 = s2(1:4);
    xyz2 = s2(5:7);
    % right of s1
    phi3 = s3(1:4);
    xyz3 = s3(5:7);

    
%     % panels of s2
%     nr = Rz(-30*pi/180)*Rz(phi2(3)*1e-9)*Ry(phi2(2)*1e-9)*Rx(phi2(1)*1e-9)*[1;0;0];
%     nl = Rz(30*pi/180)*Rz(phi2(3)*1e-9)*Ry(phi2(2)*1e-9)*Rx(phi2(1)*1e-9)*[1;0;0];
%     
%     % laser direction of s1
%     dr = Rz(30*pi/180)*Rz(phi1(3)*1e-9)*Ry(phi1(2)*1e-9)*Rx(phi1(1)*1e-9)*[1;0;0];
%     % laser direction of s3
%     dl = Rz(-30*pi/180)*Rz(phi3(3)*1e-9)*Ry(phi3(2)*1e-9)*Rx(phi3(1)*1e-9)*[1;0;0];
%     
%     offset_r = xyz1-xyz2;
%     offset_l = xyz3-xyz2;
%     
%     tr = dot(nr,offset_r)/dot(nr,dr);
%     tl = dot(nl,offset_l)/dot(nl,dl);
    
    %Jr = norm(offset_r - dr*tr)^2 - dot(nr,dr)^2;
    %Jl = norm(offset_l - dl*tl)^2 - dot(nl,dl)^2;
    siga = sig;%inv(50*eye(3));
    %sigp = 10;
    
    x21 = phi2-phi1;
    r21 = -(dphi2-dphi1);
    x23 = phi2-phi3;
    r23 = -(dphi2-dphi3);
%     if (x21-r21)<1e-2
%         factor = 1;
%     else
%         factor = 1;
%     end
    Jr = 0.5*((x21-r21))'*siga*((x21-r21))+0*exp(0.1/(155)^2*( (x21-r21)'*(x21-r21) -152^2) );
    Jl = 0.5*((x23-r23))'*siga*((x23-r23))+0*exp(0.1/(155)^2*( (x23-r23)'*(x23-r23) -152^2) );
    
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