function J = Cost3(s1,s2,s3)
    % right of s2
    phi1 = s1(1:3)*1e-9;
    xyz1 = s1(4:6);
    % s2
    phi2 = s2(1:3)*1e-9;
    xyz2 = s2(4:6);
    % s3
    phi3 = s3(1:3)*1e-9;
    xyz3 = s3(4:6);
    
    % panels of s2
    nr = Rz(-30*pi/180)*Rz(phi3(3))*Ry(phi3(2))*Rx(phi3(1))*[1;0;0];
    nl = Rz(30*pi/180)*Rz(phi3(3))*Ry(phi3(2))*Rx(phi3(1))*[1;0;0];
    
    % laser direction of s1
    dr = Rz(30*pi/180)*Rz(phi2(3))*Ry(phi2(2))*Rx(phi2(1))*[1;0;0];
    % laser direction of s3
    dl = Rz(-30*pi/180)*Rz(phi1(3))*Ry(phi1(2))*Rx(phi1(1))*[1;0;0];
    
    offset_r = xyz2-xyz3;
    offset_l = xyz1-xyz3;
    
    tr = dot(nr,offset_r)/dot(nr,dr);
    tl = dot(nl,offset_l)/dot(nl,dl);
    
    %Jr = norm(offset_r - dr*tr)^2 - dot(nr,dr)^2;
    %Jl = norm(offset_l - dl*tl)^2 - dot(nl,dl)^2;
   
    siga = 150*10^-6*eye(3);
    sigp = 10;
    c32 = [0,0,(150-30)*pi/180]';
    c31 = [0,0,(-90-30)*pi/180]';
    %disp(1*exp(-0.5*(c31-(phi1-phi3))'*sig^-1*(c31-(phi1-phi3))))
    Jr = 0*exp(-0.5*(offset_r - dr*tr)'*sigp^-1*(offset_r - dr*tr)) + 1*exp(-0.5*(c32-(phi2-phi3))'*siga^-1*(c32-(phi2-phi3)));
    Jl = 0*exp(-0.5*(offset_l - dl*tl)'*sigp^-1*(offset_l - dl*tl)) + 1*exp(-0.5*(c31-(phi1-phi3))'*siga^-1*(c31-(phi1-phi3)));
    
    l = [100,0,0]';
    s3 = Rz(210*pi/180)*l;
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