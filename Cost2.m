function J = Cost2(s1,s2,s3)
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
    nr = -Rz(-30*pi/180)*Rz(phi2(3))*Ry(phi2(2))*Rx(phi2(1))*[1;0;0];
    nl = -Rz(30*pi/180)*Rz(phi2(3))*Ry(phi2(2))*Rx(phi2(1))*[1;0;0];
    
    % laser direction of s1
    dr = Rz(30*pi/180)*Rz(phi1(3))*Ry(phi1(2))*Rx(phi1(1))*[1;0;0];
    % laser direction of s3
    dl = Rz(-30*pi/180)*Rz(phi3(3))*Ry(phi3(2))*Rx(phi3(1))*[1;0;0];
    
    offset_r = xyz2-xyz1;
    offset_l = xyz2-xyz3;
    
    tr = dot(nr,offset_r)/dot(nr,dr);
    tl = dot(nl,offset_l)/dot(nl,dl);
    
    Jr = norm(offset_r - dr*tr)^2 - dot(nr,dr)^2;
    Jl = norm(offset_l - dl*tl)^2 - dot(nl,dl)^2;

    J = Jr + Jl;
end