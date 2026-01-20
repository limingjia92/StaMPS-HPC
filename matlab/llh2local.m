function xy=llh2local(llh,origin)
% LLH2LOCAL 
%
%   Converts from longitude and latitude to local coorindates
%   given an origin.  llh (lon; lat; height) and origin should
%   be in decimal degrees. Note that heights are ignored and
%   that xy is in km.

%-------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Sept 7, 2000  Peter Cervelli        Original Code
%   Oct 20, 2000  Jessica Murray        Changed name from DM_llh2local to 
%                                       llh2local for use with non-DM functions.
%   Dec. 6, 2000  Jessica Murray        Clarified help to show that llh 
%                                       is a column vector
%   Mar 27, 2014  Andy Hooper           Correct bug for case latitude=0 
%   Jan 15, 2026  Mingjia Li            Fixed Ellipsoid Constant Bug and Added Input Dimension Check
%                                          Original code used the numeric value of 
%                                          2nd eccentricity for 1st eccentricity formulas.
%                                          Now derived correctly from WGS84 flattening (f).
%-------------------------------------------------------------

% --- Input Dimension Check  ---
% Ensure llh is 2xN or 3xN. If Nx2 or Nx3 is passed, transpose it.
if size(llh, 1) ~= 2 && size(llh, 1) ~= 3
    if size(llh, 2) == 2 || size(llh, 2) == 3
        llh = llh'; 
    else
        error('Error: Input ''llh'' must be dimensions 2xN or 3xN.');
    end
end

% --- Set Ellipsoid Constants ---
% Based on standard WGS84 definition
a = 6378137.0;              % Semi-major axis
f = 1 / 298.257223563;      % Flattening
e2 = 2*f - f^2;             % First eccentricity squared (e^2)

% --- Convert to Radians ---
llh = double(llh) * (pi/180);
origin = double(origin) * (pi/180);

% --- Do the Projection ---

% Filter indices where latitude is not zero (to avoid singularity)
z = abs(llh(2,:)) > 1e-8; % Use a small threshold instead of strict 0 for float safety

dlambda = llh(1,z) - origin(1);

% Calculate Meridian Arc Length (M) 
M = a * ((1 - e2/4 - 3*e2^2/64 - 5*e2^3/256) * llh(2,z) - ...
         (3*e2/8 + 3*e2^2/32 + 45*e2^3/1024) * sin(2*llh(2,z)) + ...
         (15*e2^2/256 + 45*e2^3/1024) * sin(4*llh(2,z)) - ...
         (35*e2^3/3072) * sin(6*llh(2,z)));

% Calculate Meridian Arc Length for Origin (M0)
M0 = a * ((1 - e2/4 - 3*e2^2/64 - 5*e2^3/256) * origin(2) - ...
          (3*e2/8 + 3*e2^2/32 + 45*e2^3/1024) * sin(2*origin(2)) + ...
          (15*e2^2/256 + 45*e2^3/1024) * sin(4*origin(2)) - ...
          (35*e2^3/3072) * sin(6*origin(2)));

% Calculate Radius of Curvature in Prime Vertical (N)
N = a ./ sqrt(1 - e2 * sin(llh(2,z)).^2);

% Coordinate calculations
E = dlambda .* sin(llh(2,z));

xy(1,z) = N .* cot(llh(2,z)) .* sin(E);
xy(2,z) = M - M0 + N .* cot(llh(2,z)) .* (1 - cos(E));

% --- Handle Special Case (Latitude near 0) ---
dlambda_eq = llh(1,~z) - origin(1);
xy(1,~z) = a * dlambda_eq;
xy(2,~z) = -M0;

% --- Convert to km ---
xy = xy / 1000;

end