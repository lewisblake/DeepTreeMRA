function R = evaluate_covariance( locs1, locs2, theta, domainGeometry )
%% COVARIANCE is a generic Covariance function
%   
%	Input: locs1, locs2, theta
% 	covariance takes two vectors of locations, and parameter vector theta
%	and calculates the the covariance function of the two
%
%	Output: R, a covariance matrix

smoothness_nu = theta(3);
% compute distances based on geometry of the data
switch domainGeometry
    case 'sphere'
    earth_radius = 6378.137; % in kilometers
    % convert from (lon, lat) to (x, y, z).
    [x1, y1, z1] = sph2cart(locs1(:,1), locs1(:,2), earth_radius);
    [x2, y2, z2] = sph2cart(locs2(:,1), locs2(:,2), earth_radius);
    % pairwise distance using chordal distance
    dist_scaled = pdist2([x1, y1, z1], [x2, y2, z2])/theta(2);
    case 'plane'
    dist_scaled = pdist2(locs1, locs2)/theta(2);   
end

if smoothness_nu == 0.5
    R = theta(1) * exp(-dist_scaled); % nu=1/2 exponential
elseif smoothness_nu == 1.5
    R = theta(1) * (1 + (sqrt(3)*dist_scaled)) * exp(-sqrt(3)*dist_scaled);
elseif smoothness_nu == 2.5
    R = theta(1) * (1 + (sqrt(5)*dist_scaled) + ((5/3)*dist_scaled^2)) * exp(-sqrt(5)*dist_scaled);
else
    not_zero = (dist_scaled ~= 0.0);
    R = theta(1) * ones(size(dist_scaled));
    R(not_zero) = theta(1) * (2^(1-smoothness_nu)/gamma(smoothness_nu)) * (dist_scaled(not_zero)).^smoothness_nu .* besselk(smoothness_nu, dist_scaled(not_zero));
end

end
