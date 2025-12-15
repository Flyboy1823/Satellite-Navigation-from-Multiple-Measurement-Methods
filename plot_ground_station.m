function Positions = plot_ground_station(GS, Plot)
%--------------------------------------------------------------------------
% plot_ground_station: plots a ground station on a spherical Earth
%--------------------------------------------------------------------------
%
% Inputs:
%   GS          - List of Nx3 of geodetic lat [deg], long [deg], alt [km]
%   markerSize  - Size of the dot (optional, default 100)
%   markerColor - Color of the dot (optional, default 'r')

lat_list = GS(:,1);
lon_list = GS(:,2);
alt_list = GS(:,3);

markerSize = 100;  % default size

% Earth's radius (mean) in km
Re = 6378;  
if Plot
    hold on
end
% Convert geodetic to ECEF (simple spherical Earth) and Plot GS
N = length(GS(:,1));
legendNames = cell(1, N);

Positions = [];

for k = 1:N
    lat = deg2rad(lat_list(k));
    lon = deg2rad(lon_list(k));
    r = Re + alt_list(k);
    x = r * cos(lat) * cos(lon);
    y = r * cos(lat) * sin(lon);
    z = r * sin(lat);
    if Plot
        legendNames{k} = sprintf('Ground Station %d', k);
        handles(k) = scatter3(x, y, z, markerSize);
    end
    Positions = [Positions; x*1000,y*1000,z*1000];
end

if Plot
    legend(handles, legendNames);
    hold off
end

end