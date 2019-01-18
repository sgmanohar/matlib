function lightbox_scatter(x,y,I, varargin)
% lightbox_scatter(x,y, I, varargin)
% x and y = vector of x/y coordinates for scatterplot -- length = size(I,4)
% I is the 4-dimensional volumes of the lesion maps, (x,y,z, subject)
% options: 
%  dim: dimension on which to flatten lesions
%  scale: size of images
%  inverty: invert the y-axis
%  brain: use brain-only template, or whole image?
% sgm 2014-2017


% 1. parse arguments
[   SCALE, INVERTY  , dim ,  BRAIN]=parsepvpairs(...
  {'scale','inverty', 'dim', 'brain' },...
  {1,      false,      1   ,  true  }, varargin{:});

% 2. load the template images:
if regexp(system_dependent('getos'), 'Linux') % where to find images?
  p1='~/Documents/MyMRIs';
  fsl_templates = '/usr/share/fsl/data/standard';
else % home windows
  p1 = '/MyMRIs/';
  if ~exist(p1,'file'), 
    p1 = 'c:\Users\smanohar\Documents\MyMRIs'; % work windows
  end
end
p1=fsl_templates;

if BRAIN
  templ = load_nii([p1 '/MNI152_T1_2mm_brain.nii.gz']);
else
  templ = load_nii([p1 '/MNI152_T1_2mm.nii.gz']);
end
templ = templ.img;

% 3. estimate scaling 
scatter(x,y); % test scatter plot to set up axes
scxlim=xlim(); scylim=ylim(); % get auto scales for normal scatter points
scale=SCALE*0.5/1000; % calculate size of voxels
% 4. plot brain images on scatter plot
for i=1:size(x,1) % for each subject
  y1=x(i);   y2=y(i); 
  if INVERTY, y2=-y2; end % invert y coords?
  % calculate image colours
  tmp = min( bsxfun( @times, sq(nansum( I(:,:,:,i), dim))', permute([1 0.05 0],[1 3 2]) ) ...
           + bsxfun( @times, sq(nanmean(templ/10000 ,dim))', permute([1 1 1],[1 3 2]) ) ...
       , 1 ) ;
  tmp2=[1 2 3]; tmp2(dim)=[]; 
  % draw image on ths scatter plot, in the appropriate location
  h=image( [-0.5, +0.5]*size(templ,tmp2(1))*scale*diff(scxlim) + y1, ...
           [-0.5, +0.5]*size(templ,tmp2(2))*scale*diff(scylim) + y2, ...
           tmp);
  set(h, 'AlphaData', mean(tmp,3)*0.8); % make outside transparent
  hold on;
end

% 5. set scaling of axes
if INVERTY, scylim= fliplr(-scylim); end
xlim(scxlim+[-50 50]*scale*diff(scxlim) ); 
ylim(scylim+[-50 50]*scale*diff(scylim) );
set(gca, 'ydir','normal')
if prod(ylim)<0    % add 'zero' lines?
  plot([0 0], ylim(),'w:'); 
end
if prod(ylim)<0
  plot(xlim(), [0 0],'w:');
end
hold off

