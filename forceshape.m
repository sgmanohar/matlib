function X = forceshape(X,S)
% force array X's shape to size S, by padding with NAN or truncating the
% array as required. It cuts of the lower or right-hand edges of the array.
SX = size(X);
baseIx = repmat({':'}, length(SX),1); % colons for each subscript
if length(SX)<length(S), 
  SX=[SX ones(1,length(S)-length(SX))];
elseif length(S)<length(SX)
  S=[S ones(1,length(SX)-length(S))];
end
for i=1:length(SX)
  SX(i)=size(X,i); % ensure the current size is correct...
  if SX(i)>S(i) % too big: truncate
    newix=baseIx;
    newix(i)={1:S(i)}; % new indices: on dimension i, take only first S(i) elements
    X=X(newix{:});
  elseif SX(i)<S(i) % too small: pad
    newix=baseIx;
    padval = nan(S(i)-SX(i),1);
    if i>1, padval = permute( padval, [2:i 1] ); end
    X=nancat(i, X, padval);
  end
end
