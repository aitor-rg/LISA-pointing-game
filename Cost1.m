function J = Cost1(s1,s2,s3)
    % s1
    phi1 = s1(1:3)*1e-9;
    xyz1 = s1(4:6);
    % left of s1
    phi2 = s2(1:3)*1e-9;
    xyz2 = s2(4:6);
    % right of s1
    phi3 = s3(1:3)*1e-9;
    xyz3 = s3(4:6);
   
    % panels of s1
    
    nr = Rz(30*pi/180)*Rz(phi1(3))*Ry(phi1(2))*Rx(phi1(1))*[1;0;0];
    
    nl = Rz(-30*pi/180)*Rz(phi1(3))*Ry(phi1(2))*Rx(phi1(1))*[1;0;0];
    
    %laser direction of s2
    dr = Rz(-30*pi/180)*Rz(phi2(3))*Ry(phi2(2))*Rx(phi2(1))*[1;0;0];
    %laser direction of s3
    dl = Rz(30*pi/180)*Rz(phi3(3))*Ry(phi3(2))*Rx(phi3(1))*[1;0;0];
    offset_r = xyz2-xyz1;
    offset_l = xyz3-xyz1;
    tr = dot(nr,offset_r)/dot(nr,dr);
    tl = dot(nl,offset_l)/dot(nl,dl);
    

    siga = 150*10^-6*eye(3);
    sigp = 10;

    %Jr = 1*exp(-0.5*(offset_r - dr*tr)'*sig^-1*(offset_r - dr*tr))-1*(dot(nr,dr)+1)^2;
    %Jl = 1*exp(-0.5*(offset_l - dl*tl)'*sig^-1*(offset_l - dl*tl))-1*(dot(nl,dl)+1)^2;
    
    c12 = [0,0,(150-(-90))*pi/180]';
    c13 = [0,0,(30-(-90))*pi/180]';
    Jr = 0*exp(-0.5*(offset_r - dr*tr)'*sigp^-1*(offset_r - dr*tr)) + 1*exp(-0.5*(c12-(phi2-phi1))'*siga^-1*(c12-(phi2-phi1)));
    Jl = 0*exp(-0.5*(offset_l - dl*tl)'*sigp^-1*(offset_l - dl*tl)) + 1*exp(-0.5*(c13-(phi3-phi1))'*siga^-1*(c13-(phi3-phi1)));

    l = [100,0,0]';
    s1 = Rz(90*pi/180)*l;
    J = Jr + Jl ;
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