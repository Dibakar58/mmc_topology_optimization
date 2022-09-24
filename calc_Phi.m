% Topology description function and derivatives
function [allPhi, allPhidrv, xval, actComp, actDsvb] = calc_Phi (allPhi, allPhidrv, xval, i, LSgrid, p, nEhcp, epsilon, actComp, actDsvb, minSz)
di = xval((i-1)*nEhcp+1: i*nEhcp) ;
x0 = di(1) ; y0 = di(2) ; z0 = di(3); l1 = di(4) + eps; l2 = di(5) + eps; l3 = di(6) + eps;
sa = sin(di(7)); sb = sin(di(8)); sg = sin(di(9));
ca = cos(di(7)); cb = cos(di(8)); cg = cos(di(9));
R = [cb*cg cb*sg -sb; sa*sb*cg-ca*sg sa*sb*sg+ca*cg sa*cb;ca*sb*cg+sa*sg ca*sb*sg-sa*cg ca*cb];
xyzLc = [LSgrid.x(:)-x0 LSgrid.y(:)-y0 LSgrid.z(:)-z0];          % local coordinates of all nodes
xyz = xyzLc*R' ;
x1 = xyz(:,1) + eps; y1 = xyz(:,2) + eps; z1 = xyz(:,3) + eps;
temp = (x1/l1).^p + (y1/l2).^p + (z1/l3).^p;
allPhi(:,i) = 1 - temp.^(1/p);                                    % TDF of i-th component
if min(di(4:6))/minSz > 0.1 && min(abs(allPhi(:,i))) < epsilon
    Ra = [0 0 0; R(3,1) R(3,2) R(3,3); -R(2,1) -R(2,2) -R(2,3)];
    Rb = [-sb*cg -sb*sg -cb; sa*cb*cg sa*cb*sg -sa*sb; ca*cb*cg ca*cb*sg -ca*sb];
    Rg = [-cb*sg cb*cg 0; -sa*sb*sg-ca*cg sa*sb*cg-ca*sg 0; -ca*sb*sg+sa*cg ca*sb*cg+sa*sg 0];
    dx1 = [[-R(1,:) 0.0 0.0 0.0]+0.0*x1 xyzLc*[Ra(1,:); Rb(1,:); Rg(1,:)]'];   % variation of x'
    dy1 = [[-R(2,:) 0.0 0.0 0.0]+0.0*y1 xyzLc*[Ra(2,:); Rb(2,:); Rg(2,:)]'];   % variation of y'
    dz1 = [[-R(3,:) 0.0 0.0 0.0]+0.0*z1 xyzLc*[Ra(3,:); Rb(3,:); Rg(3,:)]'];   % variation of z'
    temp1 = -temp.^(1/p-1).*(x1/l1).^(p-1);
    temp2 = -temp.^(1/p-1).*(y1/l2).^(p-1);
    temp3 = -temp.^(1/p-1).*(z1/l3).^(p-1);
    dpdx1 = temp1/l1; dpdy1 = temp2/l2; dpdz1 = temp3/l3;
    dpdl1 = -temp1.*(x1/l1^2); dpdl2 = -temp2.*(y1/l2^2); dpdl3 = -temp3.*(z1/l3^2);
    dpdL = 0.0*dx1; dpdL(:,4:6) =  [dpdl1 dpdl2 dpdl3];
    allPhidrv(:, nEhcp*(i-1)+1:nEhcp*i) = dpdL + dpdx1.*dx1 + dpdy1.*dy1 + dpdz1.*dz1;
else                     %deleting tiny components and removing it from active sets
    disp(['The ' sprintf('%i' ,i) '-th component is too small! DELETE it!!!']);
    allPhi(:,i) = -1e3;
    xval((i-1)*nEhcp+[4;6]) = 0;
    actComp = setdiff(actComp, i);
    actDsvb = setdiff(actDsvb, nEhcp*i-nEhcp+1:nEhcp*i);
end
end