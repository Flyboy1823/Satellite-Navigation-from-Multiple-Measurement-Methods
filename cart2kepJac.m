function J = cart2kepJac(r,v)
%--------------------------------------------------------------------------
% Computes the numerical jacobian from cartisian coordiantes to orbital
% elemnts
%--------------------------------------------------------------------------
% Inputs are Variation but should be:
%   r     - 3x1 initial position [m] in cartesian
%   v     - 3x1 initial velocity [m/s] in cartesian
% Output:
%   J  - 6x6 Jacobian Matrix
%
    x = [r; v];
    kep0=zeros(1,6);
    [kep0(1), kep0(2), kep0(3), kep0(4), kep0(5), kep0(5), kep0(6)]= ijk2keplerian(r,v);
    J = zeros(6,6);
    for j = 1:6
        for k = 1:6
            eps = 1e-6*x(k);
            dx = zeros(6,1); 
            dx(k) = eps;
            kep_plus=zeros(1,6);
            [kep_plus(1), kep_plus(2), kep_plus(3), kep_plus(4), kep_plus(5), kep_plus(6)] = ijk2keplerian(x(1:3)+dx(1:3), x(4:6)+dx(4:6));
            J(j,k) = (kep_plus(j) - kep0(j))/eps;
        end
    end
end
