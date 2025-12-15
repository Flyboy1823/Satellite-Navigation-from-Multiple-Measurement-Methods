function plotOrbitVisibility(a_fused,e_fused,inc_fused,RAAN_fused,w_fused, nu_fused, gsECI)
%--------------------------------------------------------------------------
% plotOrbitVisibility
%--------------------------------------------------------------------------
% Inputs:
%   meas = initial keplerian elements [a, e, i, RAAN, w, nu]
%   groundStations = struct array:
%       groundStations(j).lat   = [deg]
%       groundStations(j).lon   = [deg]
%       groundStations(j).alt   = [km]
%   mu   = gravitational parameter
%   t0   = initial time (seconds)
%   tf   = final time (seconds)
%   dt   = propagation step (seconds)
%

% Convert initial orbital elements to ECI
k = linspace(0, 360, 720*2);
r_test = zeros(length(k), 3);
v_test = zeros(length(k), 3);

for counter = 1:length(k)
    [r(counter, :), v(counter, :)] = keplerian2ijk(a_fused, e_fused, inc_fused, RAAN_fused, w_fused, k(counter));
end

% Convert ground stations to ECEF/ECI
numGS = 2;
N=length(k);
% Determine LoS visibility for each point
vis1 = false(N,1);
vis2 = false(N,1);

for i = 1:length(k)
    [el1, ~] = elevazim(r(i,:), gsECI(1,:));
    [el2, ~] = elevazim(r(i,:), gsECI(2,:));

    vis1(i) = el1 > 0;
    vis2(i) = el2 > 0;
end

% Four visibility cases
idx_none = ~vis1 & ~vis2;
idx_gs1  =  vis1 & ~vis2;
idx_gs2  = ~vis1 &  vis2;
idx_both =  vis1 &  vis2;

% Split r into 4 sets
r_none = r(idx_none, :);
r_gs1  = r(idx_gs1,  :);
r_gs2  = r(idx_gs2,  :);
r_both = r(idx_both, :);

% Assign Colors
% None visible  → black
% GS1 only      → blue
% GS2 only      → red
% Both          → green

% colors = zeros(N,3);
% 
% for k = 1:N
%     if vis1(k) && vis2(k)
%         colors(k,:) = [0,1,0];         % green
%     elseif vis1(k)
%         colors(k,:) = [0,0,1];         % blue
%     elseif vis2(k)
%         colors(k,:) = [1,0,0];         % red
%     else
%         colors(k,:) = [0,0,0];         % black
%     end
% end

% Plot the Earth
figure
hold on
earth_sphere('km');
axis equal
hold on

hNone = scatter3(r_none(:,1)/1000, r_none(:,2)/1000, r_none(:,3)/1000, 10, 'k');
hGS1  = scatter3(r_gs1(:,1)/1000,  r_gs1(:,2)/1000,  r_gs1(:,3)/1000, 10, 'b');
hGS2  = scatter3(r_gs2(:,1)/1000,  r_gs2(:,2)/1000,  r_gs2(:,3)/1000, 10, 'r');
hBoth = scatter3(r_both(:,1)/1000, r_both(:,2)/1000, r_both(:,3)/1000, 10, [0.5 0 0.5]); % purple

hold on
earth_sphere('km');
axis equal

hold on
hgs1=scatter3(gsECI(1,1)/1000, gsECI(1,2)/1000, gsECI(1,3)/1000, 200, 'b', 'filled');  % GS1 color
hgs2=scatter3(gsECI(2,1)/1000, gsECI(2,2)/1000, gsECI(2,3)/1000, 200, 'r', 'filled');  % GS2 color



xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Orbit Visibility by Ground Stations');
legend([hNone hGS1 hGS2 hBoth hgs1 hgs2], ...
       {'No LOS','LOS GS1 Only','LOS GS2 Only','LOS Both', 'GS1 Station', 'GS2 Station'}, ...
       'Location','bestoutside');
grid on;
hold off;

end
