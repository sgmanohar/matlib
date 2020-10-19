function stat = vlsm(I, X, varargin)
% stat = vlsm(I, X, ...)
% voxel-based lesion-symptom mapping
% regress behaviour against lesion degree in each voxel.
% I = imaging ( x,y,z, subjects )
% X = design matrix ( subject, regressor )
%     each column is run in a separate VLSM.
% RANKED      = 0   % rank transform before t-test
% PERMUTATION = 0   % run permutation test instead of plain t test
% FILTER_N    = 4   % how many patients need to have lesions in a voxel to warrant analysing it
% ADD_CONTROL = []  % add healthy control values as dummy 'no-lesion' subjects (matrix with same #cols as X) 
% COVARIATES  = []  % values to regress out before VLSM (matrix with same #rows as X  
% PRESMOOTH   = 0   % num voxels gaussian kernel smoothing before vlsm
% method      = 'residual' (default) or 'vlsm' - way of correcting for covariates 
% mirror      = false  % treat L/R lesions as comparable
% savepngfile = the filename prefix. Do not save if [] or false, 
% sgm 2014

%% smoothing required

[      RANKED, PERMUTATION,   MIRROR,   FILTER_N, ADD_CONTROL, PRESMOOTH, COVARIATES, method,   SAVE_PNG  ]...
  = parsepvpairs( {
      'RANKED', 'PERMUTATION','MIRROR','FILTERN','ADDCONTROL','PRESMOOTH','COVARIATES','METHOD',  'SAVEPNGFILE' 
  },{  0       0                0       4           []           0          []         'residual', 'vlsm'      
  }, varargin{:});
% 
Plot = enum({'P','STAT', 'B'});

PLOT = Plot.P;

if isstruct(X)
  f=fieldnames(X);
  Xd=[];
  for i=1:length(f)
    Xd=[Xd X.(f{i})];
  end
  regressor_names = f;
else
  regressor_names = arrayfun(@(x)sprintf('X%g',x), 1:size(X,2),'uni',0);
end

sz = size(I);
Is = double(I); % convert from int16 to double, needed for OLS
if PRESMOOTH
  % this is the FWHM as passed to spm_smooth
  Is = smooth_volume_using_spm( double(Is), PRESMOOTH );
end

if MIRROR % add on X-mirror image to each lesion map
  % If = I flattened
  If = reshape( Is + Is(end:-1:1,:,:,:) , [], sz(4) )';
else
  If = reshape( Is, [], sz(4) )'; % subjects x 902629 voxels
end


if RANKED && ~exist('Ifr','var') % calculate rank of subjects for each voxel?
  for i=1:size(If,2) % calculate a rank-ordered Y. For each voxel
    if any(If(:,i)>0) % any people with lesions at this voxel?
      Yr(:,i) = rankorder(If(:,i)); % replace the lesion value with its rank
    else % otherwise everyone's rank is zero. 
      Yr(1:sz(4),i) = 0;
    end
  end
  Ifr=int8(Yr); % voxelwise-ranked image
end

% which regerssors to calculate?
selection  = 1:size(X,2);  
% selct areas with enough variance FILTER(VOXEL) for permutation test
filter_roi_f = sum( If>0 ) >= FILTER_N; 
for j=selection % for each regressor, 
  tic;
  % get the regressors
  regressor  = X(:,j);
  if ADD_CONTROL
    cregressor = ADD_CONTROL;         % matrix of regressors for controls
    NC         = length(cregressor);  % number of controls to add
  else
    cregressor = [];
  end
  %%%%%%% added Oct 2015: factor out covariates before VLSM
  lesvol     = sum(reshape(I,[],sz(4)))'; 
  covariate  = []; % regressor matrix, will eventually be   length(include) x length(covariates)
  for i=1:size(COVARIATES,2) % for each covariate to factor out, 
    cov_i = nanzscore(COVARIATES(:,i));
    cov_i(isnan(cov_i))=0; % remove nans from the covariate
    covariate = [covariate cov_i]; % add it to the regressor matrix
  end
  switch method
    case 'residual'
      if ~isempty(covariate) % remove covariates by regression, and use only residuals
        [~,~,regressor] = regress(regressor, [ones(size(regressor)) covariate]);
      end % now 'regressor' will be decorrelated from the specified covariates.
  end
  %%%%%%
  
  regressor  = nanzscore(regressor);   
  cregressor = nanzscore(cregressor);
  all_regressor{j}  = regressor;   % store for later use
  all_cregressor{j} = cregressor;
  
  % create design matrix
  Xd = [ ones(length(regressor),1) , regressor ];
  switch method
    case 'residual'
      contrast = [0 1];
    case 'vlsm'
      Xd = [Xd covariate]; contrast = [0 1 zeros(1,size(covariate,2))];
  end
  % exclude any subjects with no regressor value
  exclude = find(isnan(regressor)); 
  Y = If; % load Y afresh
  Xd (exclude,:)=[]; Y (exclude,:)=[]; % remove exclusions
  if ADD_CONTROL % for controls add lesion mask of zeros 
    Xd = [Xd; [ones(NC,1), cregressor ] ];
    Y = cat(1, Y, zeros([NC size(Y,2)]) );   
  end
  if RANKED % similarly prepare ranked dataexclude = 16; % exclude blackshaw
    Xr = [Xd(:,1) rankorder(Xd(:,2))']; % ranked version of X
    Yr=Ifr;
    if ADD_CONTROL
      % the controls all have 'rank zero' lesion; so they take the rank of
      % the middle one. **This might not be quite right though, since many 
      % voxels in patients also have zero lesion....
      Yr=[Yr+NC; ones([NC size(Yr,2)]) * NC/2]; 
    end
    Xr(exclude,:)=[]; Yr(exclude,:)=[];
  end

  % Least squares regression
  [b0,~,t0] = ols( Y, Xd, contrast );
  % [r0]      = corr( Y,X(:,2), 'type','spearman' );
  if RANKED
    [r0,~,t0] = ols( double(Yr), Xr, contrast );
  end
  
  if PERMUTATION % DO PERMUTATIONS ?
    NI = 5000; % how many permutations
    bmax  = nan(1,NI);   bmin = nan(1,NI);
    tmax  = nan(1,NI);   tmin = nan(1,NI);
    brmax = nan(1,NI);  brmin = nan(1,NI);
    trmax = nan(1,NI);  trmin = nan(1,NI);
    for i=1:NI  % for each permutation
      perm = randperm(size(Y,1));
      [ bi, ~, ti ] = ols( Y, Xd(perm,:), contrast );
      % [ ri ] = corr( Y, X(perm,2), 'type','spearman')
      bmax(i) = max(bi(filter_roi_f)); % build distribtion of maximum statistic over all voxels
      bmin(i) = min(bi(filter_roi_f));
      tmax(i) = max(ti(filter_roi_f));
      tmin(i) = min(ti(filter_roi_f));
      p_t0 = tcdf(t0,size(Y,1));
      if RANKED % do RANKED versions?
        [ bri, ~, tri ] = ols( Yr, Xr(perm,:), contrast );
        brmax(i) = max(bri(filter_roi_f));
        brmin(i) = min(bri(filter_roi_f));
        trmax(i) = max(tri(filter_roi_f));
        trmin(i) = min(tri(filter_roi_f));
      end
      if mod(i,100)==0; fprintf('.'); end; % could take some time, --> user feedback
    end
    
    % calculate P such that P>0.95 is significant
    % upper tail
    if ~isunix % windows version can handle bsxfun
      stat(j).p_bu = reshape( mean(bsxfun(@lt, bmax, b0'),2),  sz(1:3) ); % P( b0 > bmax )?
      stat(j).p_tu = reshape( mean(bsxfun(@lt, tmax, t0'),2),  sz(1:3) );
      % lower tail
      stat(j).p_bl = reshape( mean(bsxfun(@gt, bmin, b0'),2),  sz(1:3) ); % P( b0 < bmin )?
      stat(j).p_tl = reshape( mean(bsxfun(@gt, tmin, t0'),2),  sz(1:3) ); % proportion of permutations which have min(t) above real t.
      if RANKED
        stat(j).p_bru = reshape( mean(bsxfun(@lt, brmax, r0'),2),  sz(1:3) );
        stat(j).p_tru = reshape( mean(bsxfun(@lt, trmax, r0'),2),  sz(1:3) );
        stat(j).p_brl = reshape( mean(bsxfun(@gt, brmin, r0'),2),  sz(1:3) );
        stat(j).p_trl = reshape( mean(bsxfun(@gt, trmin, r0'),2),  sz(1:3) );
      end
    else % unix version wastes tonnes of memory when doing mean(bsxfun)
      for i=1:length(b0) % so do it manually - for each voxel
        tmp_p_bu(i) = mean( bmax<b0(i) );
        tmp_p_tu(i) = mean( tmax<t0(i) );
        tmp_p_bl(i) = mean( bmin>b0(i) );
        tmp_p_tl(i) = mean( tmin>t0(i) );
      end
      stat(j).p_bu = reshape( tmp_p_bu, sz(1:3) );
      stat(j).p_tu = reshape( tmp_p_tu, sz(1:3) );
      stat(j).p_bl = reshape( tmp_p_bl, sz(1:3) );
      stat(j).p_tl = reshape( tmp_p_tl, sz(1:3) );
    end
    lightbox_images_for_thesis( stat(j).p_bu.*(stat(j).p_bu>0.95) );
  else % NO PERMUTATION?
    % use simple t-tests uncorrectedindependent of the 4 
    p_t0 = tcdf(t0,size(Y,1)-1);
    p_t0(p_t0<0.5) =  -p_t0(p_t0<0.5);  % mirror so that values between -0.05 and +0.05 are significant
    p_t0(p_t0>0.5) = 1-p_t0(p_t0>0.5); % and sign reflects the direction of significance.
    stat(j).p_t = reshape( p_t0,  sz(1:3) );
    stat(j).t   = reshape( t0,  sz(1:3) );
    stat(j).b   = reshape( b0,  sz(1:3) );
  end
  time(j)=toc;
end

if PERMUTATION && 1
  save(sprintf('vlsm_permuted_ols_results'), 'stat', 'PRESMOOTH','FILTER_N','MIRROR','ADD_CONTROL');
end


switch PLOT
%% image of each statistic in stat?
  case Plot.STAT
    k=1; clear tmp
    ss=fieldnames(stat);
    for i=1:length(ss); % statistic
      for j=1:size(X,2)-1 % beh_var
        tmp(:,:,:,k)=stat(j).(ss{i});
        name{k} = sprintf('regressor %g\n%s',j, ss{i});
        k=k+1;
      end
    end
    lightbox_images_for_thesis(tmp.*(tmp>0.95), name, 'dim',3);
    any_signif = sq(sum(sum(sum(tmp>0.95))))

%% this seems to produce pretty pictures
% after running VLSM, we have stat(1) for velocity reward sensitivity
% with b and t statistics.
% note that this version draws all statistics, significant or not.
  case Plot.B
    b0=stat(1).t;
    tmp2=reshape( b0,  sz(1:3) );  % get the basic b for each voxel. (coefficient of vel rew sens)
    figure(1); colormap(hot(255)); % where is b positive?
    % smoothing of 2, red colours
    smth = @(x) smooth_volume_using_spm( x, 2 );
    lightbox_images_for_thesis( smth( tmp2.*(tmp2>0) ), regressor_names, 'Dim',3, 'Invert',0)
    figure(2); colormap(flipud(1-hot(255))) % where is b negative? blue colours
    lightbox_images_for_thesis( smth( -tmp2.*(tmp2<0) ), regressor_names, 'Dim',3, 'Invert',0)

%% but this version is correct i.e. only draws what is significant
% now filter based on power
% variance in each voxel:
  case Plot.P
    include  = 1:sz(4);
    powermap = var(cat(4, Is(:,:,:,include), Is(end:-1:1,:,:,include)),[],4);
    % how many patients in each voxel
    numbermap= (sum( (Is(:,:,:,include)>0.01) + (Is(end:-1:1,:,:,include)>0.01) ,4) );
    MINPOWER = 0.0;  % exclude voxels with low lesion variance
    MINPOWERN= 4;    % exclude voxels with subjects fewer than this. Should be same as FILTER_N
    ALPHA    = 0.05;  % p-value threshold (two-tailed = 0.05, one tailed=0.10)
    POSTSMOOTH   = 5;     % final image smooth fwhm
    DIM      = 1;
    if 0 % filter by POWER
      % filter by minimum variance
      powerfilter = repmat(powermap >=MINPOWER,  [1 1 1 length(selection)]);
    else % filter by NUM LESIONS
      % filter by number of patients in each voxel
      powerfilter = repmat(numbermap>=MINPOWERN, [1 1 1 length(selection)]);
    end
    figure(1);
    if PERMUTATION
      x=cat(4, stat(selection).p_tu); x=1-x; % value close to 1 is significant, so invert x on interval 0-1
    else
      x=cat(4, stat(selection).p_t); % value close to 0 is significant
    end
    x=(x<ALPHA & x>0).*(-log(abs(x)+eps)) .* powerfilter;
    significant_voxels = sq(sum(sum(sum(x>0))));
    x=smooth_volume_using_spm(x,POSTSMOOTH);
    lightbox_images_for_thesis( x , regressor_names ,'invert',0, 'dim', DIM)
    xupper=x;
    figure(2)
    if PERMUTATION
      x=cat(4, stat(selection).p_tl);  x=1-x; % value close to 0 is significant so invert x on interval 0-1
    else
      x=cat(4, stat(selection).p_t); x=-x; % negative values close to 0 are significant, so negate x
    end
    x=(x<ALPHA & x>0).*(-log(abs(x)+eps)) .* powerfilter;
    significant_voxels(:,2) = sq(sum(sum(sum(x>0))))
    x=smooth_volume_using_spm(x,POSTSMOOTH);
    lightbox_images_for_thesis( x , regressor_names ,'invert',0, 'dim', DIM)
    xlower=x;
    if SAVE_PNG
      figure(1); saveas(gcf,sprintf('%s_positive_mirror%g_smooth%g_numthresh%g_presmooth%g_perm%g_p%g.png',...
        SAVE_PNG, MIRROR, POSTSMOOTH, MINPOWERN, PRESMOOTH, PERMUTATION, ALPHA),'png');
      figure(2); saveas(gcf,sprintf('vlsm_negative_mirror%g_smooth%g_numthresh%g_presmooth%g_perm%g_p%g.png',...
        SAVE_PNG, MIRROR, POSTSMOOTH, MINPOWERN, PRESMOOTH, PERMUTATION, ALPHA),'png');
    end
end

return 



%%%%%%%%%%%%%%%%%%
% first get amount of lesion at the locus of interest, for each patient
xyz_mni = peaks_mni{peak_ix,1}; % this array is defined in AnalyseTrioPatients/2014/myvlsm / region_of_interest.m
xyz=MNI(xyz_mni);  % this function is in region_of_interest.m
test=zeros(sz(1:3));test(xyz(1),xyz(2),xyz(3))=1;
RADIUS  = 5; % (exactly as used in the main analysis)
voi_xyz=smooth_volume_using_spm( test, RADIUS );
amtlesion_xyz = sq(sum(sum(sum( bsxfun(@times, voi_xyz, Ic(:,:,:,include)) )))); 


