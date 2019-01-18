function O=rotatend(theta)
% O=rotatend(theta)
% Rotate in n-dimensions
% theta is a vector of angles, where element i corresponds to rotation in
% plane in dimensions [i,i+1]. Rotations are composed in order of ascending
% dimensions.
% sgm 2016
D = length(theta); % dimensionality of manifold
O = eye(D); % initialise with identity
for i=1:length(theta) % for each angle
  % corresponding 2-dimensional rotation
  R=[  cos(theta(i)), sin(theta(i)) zeros(1,D-2)
      -sin(theta(i)), cos(theta(i)) zeros(1,D-2)
      zeros(D-2,2)                  eye(D-2)     ];
  ord=circshift([1:D]', i-1);
  R = R(ord,ord);
  O=O*R; % compose with current transform
end
  