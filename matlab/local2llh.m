function llh = local2llh(xy, origin)
% LOCAL2LLH 
%
%   Converts from local coordinates (x,y in km) to longitude and latitude 
%   given an origin [lon, lat] in decimal degrees. Heights are ignored.
%   This is an iterative solution for the inverse of a polyconic projection.
%   Output is [lon; lat] in decimal degrees.
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%   Sep 7, 2000   Peter Cervelli        Original Code
%   Apr 4, 2001   Peter Cervelli        Added failsafe to avoid infinite loop
%   Aug 23, 2001  Jessica Murray        Clarification to help.
%   Feb 2026      Mingjia Li            Fixed Ellipsoid Constant Bug, optimized
%                                       code, and added dimension check.
%                                       (Original used 2nd eccentricity numeric
%                                       value for 1st eccentricity formulas).
%   ======================================================================

% --- Input Dimension Check ---
% Ensure xy is 2xN. If Nx2 is passed, transpose it.
if size(xy, 1) ~= 2
    if size(xy, 2) == 2
        xy = xy'; 
    else
        error('Error: Input ''xy'' must be dimensions 2xN.');
    end
end

% --- Set Ellipsoid Constants (WGS84) ---
a = 6378137.0;              % Semi-major axis
f = 1 / 298.257223563;      % Flattening
e2 = 2*f - f^2;             % First eccentricity squared (e^2)

% --- Convert units (km to meters, degrees to radians) ---
xy = double(xy) * 1000;
origin = double(origin) * (pi/180);

% --- Perform Inverse Projection ---
% Calculate Meridian Arc Length for Origin (M0)
M0 = a * ((1 - e2/4 - 3*e2^2/64 - 5*e2^3/256) * origin(2) - ...
          (3*e2/8 + 3*e2^2/32 + 45*e2^3/1024) * sin(2*origin(2)) + ...
          (15*e2^2/256 + 45*e2^3/1024) * sin(4*origin(2)) - ...
          (35*e2^3/3072) * sin(6*origin(2)));

% Filter indices where latitude is not exactly zero (Float safe)
z = abs(xy(2,:) + M0) > 1e-8;

% Initialize llh array (2xN)
llh = zeros(2, size(xy, 2));

% Pre-calculate constants for iteration
A = (M0 + xy(2,z)) / a;
B = xy(1,z).^2 ./ a^2 + A.^2;

llh(2,z) = A;
delta = Inf;
c = 0;

% --- Iterative Solution Loop ---
while max(abs(delta)) > 1e-8
    
    C = sqrt(1 - e2 * sin(llh(2,z)).^2) .* tan(llh(2,z));

    M = a * ((1 - e2/4 - 3*e2^2/64 - 5*e2^3/256) * llh(2,z) - ...
             (3*e2/8 + 3*e2^2/32 + 45*e2^3/1024) * sin(2*llh(2,z)) + ...
             (15*e2^2/256 + 45*e2^3/1024) * sin(4*llh(2,z)) - ...
             (35*e2^3/3072) * sin(6*llh(2,z)));

    Mn = 1 - e2/4 - 3*e2^2/64 - 5*e2^3/256 ...
         - 2 * (3*e2/8 + 3*e2^2/32 + 45*e2^3/1024) * cos(2*llh(2,z)) ...
         + 4 * (15*e2^2/256 + 45*e2^3/1024) * cos(4*llh(2,z)) ...
         - 6 * (35*e2^3/3072) * cos(6*llh(2,z));

    Ma = M / a;

    % Newton-Raphson Delta
    delta = -(A .* (C .* Ma + 1) - Ma - 0.5 * (Ma.^2 + B) .* C) ./ ...
            (e2 * sin(2*llh(2,z)) .* (Ma.^2 + B - 2 * A .* Ma) ./ (4 * C) + ...
            (A - Ma) .* (C .* Mn - 2 ./ sin(2*llh(2,z))) - Mn);

    llh(2,z) = llh(2,z) + delta;

    c = c + 1;
    if c > 100
        error('Convergence failure in local2llh inverse projection.');
    end
end

% Finalize Longitude calculation
llh(1,z) = (asin(xy(1,z) .* C / a)) ./ sin(llh(2,z)) + origin(1);

% --- Handle Special Case (Latitude near 0) ---
llh(1,~z) = xy(1,~z) / a + origin(1);
llh(2,~z) = 0;

% --- Convert back to decimal degrees ---
llh = llh * (180/pi);

end