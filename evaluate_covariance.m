function R = evaluate_covariance( locs1, locs2, theta )
%% COVARIANCE is a generic Covariance function
%   
%	Input: locs1, locs2, theta
% 	covariance takes two vectors of locations, and parameter vector theta
%	and calculates the the covariance function of the two
%
%	Output: R, a covariance matrix

% Calculate the pairwise distance, h, and divide (scale) by the second entry of theta
% Let the covariance matrix be the first parameter of theta mutliplied by an exponential
% determined by h
R = theta(1) * exp(-pdist2(locs1, locs2)/theta(2)); % nu=1/2 exponential

end