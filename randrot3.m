function M=randrot3


W = randn(4,1);     % random quaternion
W = W/norm(W);      % make length = 1
W = num2cell(W);    % convert to w,x,y,z
[w,x,y,z] = deal(W{:});
M = [
  w 	z 	-y 	x
  -z 	w 	x 	y
  y 	-x 	w 	z
  -x 	-y 	-z 	w
]	* [
  w 	z 	-y 	-x
  -z 	w 	x 	-y
  y 	-x 	w 	-z
  x 	y 	z 	w
];
