function generalRaceFitExample()


% example data
test=@(N)1./(3+randn(N,1));
%%%%   RESP1_T              RESP2_T
t = [ test(80)         inf*ones(80,1)   % condition 1 correct
      inf*ones(20,1)   2*test(20)       % condition 1 error
      test(80)         inf*ones(80,1)   % condition 2 correct
      inf*ones(20,1)   2*test(20)       % condition 2 error
    ];

% show example data
reciprobit(t,  [ repmat([1 2],100,1) 
                 repmat([3 4],100,1) ] ,  [0.80 0.20 0.80 0.20] )


% set up a race model
raceParams(1).meanRate  = [1 1];
raceParams(2).meanRate  = [1 1];
raceParams(1).condition = 1;
raceParams(2).condition = 2;

% do simulated annealing!
generalRaceFit2(t, [ 1*ones(100,1)
                     2*ones(100,1) ],... condition 1 or 2
                raceParams, @varyParamFunc, [1 1 1]);
              
              
% set parameters:
%  try to fit the common stdRate, and 
%  independently fit the two conditions' correct meanRate of rise.
function p=varyParamFunc(p,x)
    p.stdRate = x(1);
    if     p.condition == 1,  p.meanRate(1) = x(2);
    elseif p.condition == 2,  p.meanRate(1) = x(3);
    end