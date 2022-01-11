function b_rot = stewartrot(theta, phi, psi, b_home)
    rotmat = zeros(4) ;
    rotmat(1,1) = cosd(psi)*cosd(phi) ;
    rotmat(1,2) = cosd(psi)*sind(phi)*sind(theta) - sind(psi)*cosd(theta) ;
    rotmat(1,3) = cosd(psi)*sind(phi)*cosd(theta) + sind(psi)*sind(theta) ;
    rotmat(1,4) = 0 ;
    rotmat(2,1) = sind(psi)*cosd(phi) ;
    rotmat(2,2) = sind(psi)*sind(phi)*sind(theta) + cosd(psi)*cosd(theta) ;
    rotmat(2,3) = sind(psi)*sind(phi)*cosd(theta) - cosd(psi)*sind(theta) ;
    rotmat(2,4) = 0 ;
    rotmat(3,1) = -sind(phi) ;
    rotmat(3,2) = cosd(phi)*sind(theta) ;
    rotmat(3,3) = cosd(phi)*cosd(theta) ;
    rotmat(3,4) = 0 ;
    rotmat(4,:) = [0, 0, 0, 1] ;
    b_rot = rotmat*b_home ;
end
