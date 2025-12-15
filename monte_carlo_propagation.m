function [x_mc, P_mc] = monte_carlo_propagation(r0, v0, P0, tprop, N)
%--------------------------------------------------------------------------
% MONTE_CARLO_PROPAGATION Propagate orbit uncertainty using Monte Carlo
%--------------------------------------------------------------------------
%
% Inputs:
%   x0    : 6x1 initial state [r0; v0]  [m, m/s]
%   P0    : 6x6 initial covariance  [m^2, (m/s)^2]
%   tprop : propagation time [s]
%   mu    : gravitational parameter [m^3/s^2]
%   N     : number of Monte Carlo particles
%
% Outputs:
%   x_mc  : propagated mean state
%   P_mc  : propagated covariance

    % Generate Monte Carlo particles
    P0 = (P0 + P0')/2;  % force symmetry
    [eigvec, eigval] = eig(P0);
    eigval(eigval < 0) = 0;  % set negative eigenvalues to zero
    P_psd = eigvec * eigval * eigvec';
    X0 = mvnrnd([r0;v0], P_psd, N)';   % 6 x N

    % Preallocate
    Xprop = zeros(size(X0));

    % Propagate each particle
    for i = 1:N
        [r_mc(:,i) v_mc(:,i), ~, ~, ~] = Propagate(tprop, X0(1:3, i), X0(4:6, i), 0);
    end

    % Compute mean and covariance
    x_mc = mean([r_mc; v_mc], 2);
    P_mc = cov([r_mc; v_mc]');
end
