function lightbox_images_for_thesis(I, pf_names, varargin)
% function lightbox_images_for_thesis(I, pf_names, ...)
% 
% 'Invert'           white background
% 'ShowNames' 
% 'AmplifyLesions'   rescale maps to [0:1] before drawing
% 'Dim'              direction to slice - 1,2,3 
% 'ChopNegative'     truncate masks to [0:1]
% 'UseColormap'      current colormap is used, otherwise uses single color
% 'TemplateColor'    color of brain image behind - [.15 .15 .15]
% 'Transparency'     can you see slices behind other slices?
% 'Transform'        after slicing the input I, pass it through a function 
%                    to transform (default = identity = @(x)x )
% 'WholeFigure'      clear figure, draw whole page view (1). If 0, then
%                    draw in current axis
% 'Range'            slice range eg [ 26:4:75 ]
% 'OneHemisphere'    show only one hemisphere (1 or 2 for left/right)
% 'Colorbar'         draw a colour bar for the values of the input image
% Create lightbox views for thesis


matlib_folder = [getExperimentBase '/matlib'];
if regexp(system_dependent('getos'), 'Linux')
  pfcdata='~/Documents/Experiment/Sanjay/Battery/AnalyseTrioPatients/2013-07-23/pfc_temp.mat';
  p1='~/Documents/MyMRIs';
else
  pfcdata='S:\Experiment\Sanjay\Battery\AnalyseTrioPatients\2013-07-23\pfc_temp.mat'; % home windows
  p1 = '/MyMRIs/';
  if ~exist(p1,'file'), 
    p1 = 'c:\Users\smanohar\Documents\MyMRIs'; % work windows
  end
  if ~exist(p1,'file')
      p1 = [ getExperimentBase() '/matlib' ];
  end
end
p1 = matlib_folder; 

fsl_templates = '/usr/share/fsl/data/standard';

if ~exist('I','var')
  % includes I
  load(pfcdata);
  sub=1:20;
  %sub(find(strcmpi(pf_names, 'blackshaw')))=[];
  % select subjects:
  sub=[1 11 3 5 2 7 6 8 4 9 10  12 13 14 17 15  18 19 20];
else
  sub=1:size(I,4);
end

% templates:
BRAIN = false;
if BRAIN
  templ = load_nii([p1 '/MNI152_T1_2mm_brain.nii.gz']);
else
  templ = load_nii([p1 '/MNI152_T1_2mm.nii.gz']);
end
templ = templ.img;
if size(templ,1) ~= size(I,1) || size(templ,2)~=size(I,2)
    [yy,xx,zz]=meshgrid( linspace(0,1,size(I,2)), ...
                         linspace(0,1,size(I,1)), ...
                         linspace(0,1,size(I,3)) );
    templ = interpn( linspace(0,1,size(templ,1)), ...
                     linspace(0,1,size(templ,2)), ...
                     linspace(0,1,size(templ,3)), double(templ),...
                     xx,yy,zz...
             );
end
%%

%sub=sub(1:9);
%sub=sub(10:19);

% contrast scaling of template
% midpoint and width of template color scale
% C = 10000 * ( C_BASE + C_SCALE / (1+e^(-(X-C_MID)/C_WID)) )
INVERT          = 1;
SHOW_NAMES      = 1 && exist('pf_names','var');
AMPLIFY_LESIONS = 1;
DIM             = 3;
CHOP_NEGATIVE   = 1; %truncate negative lesion values to zero
USE_COLORMAP    = 1;
ONE_HEMISPHERE  = 0; % show only one hemisphere?
COLORBAR        = 0; % show a color bar specifying the values of the image
WIDTH           = 800; % pixels, total width of full image, all slices.

[INVERT SHOW_NAMES AMPLIFY_LESIONS DIM CHOP_NEGATIVE USE_COLORMAP ...
  GLOBAL_OVERLAP, RESCALE_INDIVIDUALLY, RANGE, col0, transparency, ...
  transform, WHOLE_FIGURE, ONE_HEMISPHERE, COLORBAR, WIDTH ...
  DO_RESCALE] ...
  = parsepvpairs({'Invert','ShowNames','AmplifyLesions','Dim','ChopNegative','UseColormap', ...
  'GlobalOverlap', 'RescaleIndividually', 'Range', 'TemplateColor', 'Transparency' ...
  'Transform', 'WholeFigure' , 'OneHemisphere', 'Colorbar', 'Width', 'doRescale'} ...
               , { 0     ,  1        ,  1             , 3   , 1            , 1           ...
  , 0,             1 ,                  [],         [1 1 1]*0.15 , 0 ...
  @(x)x      , 1,       0,   0,   800, true} ...
               , varargin{:} );  
SHOW_NAMES      = 1 && exist('pf_names','var') && ~isempty(pf_names);


if ~INVERT
 c_mid=8000;c_wid=1300;c_scale=5; c_base=0; % noniverted contrast
 col1 = [1 0.4, 0];  % colour of mask
else
  c_mid=8000;c_wid=2000;c_scale=4.5; c_base=-0.1; % inverted contrast
  col1 = [0.2 0.7, 1];  % colour of mask
end

  

% slices
if      DIM==3, % for dim = 3: AXIAL [26:6:75]  % AXIAL
  slices = [26:4:75];  if ~isempty(RANGE), slices=RANGE; end
  slicer = {':',':',slices};
  ht = size(I,2);
elseif  DIM==2, % for dim = 2: [50:6:100]      % CORONAL
  slices = [50:4:100]; if ~isempty(RANGE), slices=RANGE; end
  slicer = {':',slices,10:91};
  ht = length(slicer{3});
elseif  DIM==1, % for dim = 1:  [31:4:60]      % SAGITTAL
  slices = [31:2:60]; if ~isempty(RANGE), slices=RANGE; end
  slicer = {slices,':',':'}; 
  slicer = {slices,[109:-1:1],':'};  % flip each slice LR so you can see frontal lesions
  ht = 91;
end

midx = floor(size(I,1)/2);
if ONE_HEMISPHERE==1 % if only showing one hemisphere, adjust image and slicing
  I     = I(1:midx+1, :,:,:);
  templ = templ(1:midx+1,:,:,:);
  % slicer{1}(slicer{1}>midx+1)=[]; % remove slices after midpoint
  if DIM==1
    slices(slices>midx+1)=[];
  end
elseif ONE_HEMISPHERE==2
  I     = I(midx:end, :,:,:);
  templ = templ(midx:end,:,:,:);
  % slicer{1}(slicer{1}<midx)=[];
  if DIM==1
    slices(slices<midx)=[];
  end
end

DX=60; % distance between slices
% WIDTH = 800;
DX = floor((WIDTH-110)/length(slices));
if exist('pf_names','var')
  pf_names=deCamel(pf_names); % make names look more friendly
end

allsubj=[];
for i=1:length(sub) % for each subject
  Image = I(:,:,:,sub(i));
  Image(isnan(Image))=0; % remove nans - they cause a whiteout mask at present.
  if GLOBAL_OVERLAP % GLOBAL OVERLAP?
    Image = sum(I(:,:,:,sub),4);
  end
  if ~AMPLIFY_LESIONS && any(Image(:)>0) && RESCALE_INDIVIDUALLY % rescale each volume
    Image=Image/max(Image(:));
  elseif ~AMPLIFY_LESIONS && any(I(:))>0 && DO_RESCALE % rescale across whole set of volumes
    Image=Image/max(I(:));
  end
  final=zeros([ht,WIDTH,3]); % overlaid output
  for j=1:length(slices)
    subplot(length(sub),1,i)
    sl=slices(j); slicer{DIM}=sl;
    les  = sq(nansum(Image( slicer{:} )  , 4))';
    les  = transform(les);
    if AMPLIFY_LESIONS
      maxlesion = max(les(:));
      if maxlesion>0 , les  = les/maxlesion; end;
    end
    if any(les(:)<0) && CHOP_NEGATIVE, les(les<0) = 0; end
    tmpl = double(sq(templ( slicer{:} )))';
    tmpl = 10000*(c_base+c_scale*1./(1+exp(-(tmpl-c_mid)/c_wid)));
    if USE_COLORMAP
      lesind = floor(max(min(les,1),0)*(size(colormap,1))+1);
      tmp  = ind2rgb( lesind,   [0 0 0;colormap] ) + ...
             bsxfun(@times, tmpl,  permute(col0, [1 3 2])/8192 );
    else
      tmp  = bsxfun(@times, les,   permute(col1, [1 3 2])/1 ) + ...
             bsxfun(@times, tmpl,  permute(col0, [1 3 2])/8192 );
    end
    tmp=max(0, min(1, tmp));
    offset=(j-1)*DX;
    % merge
    region = {1:size(tmp,1), offset+1:offset+size(tmp,2) ,1:3};
    %final( 1:size(tmp,1), offset+1:offset+size(tmp,2) ,: ) = ...
    %   final( 1:size(tmp,1), offset+1:offset+size(tmp,2),:  ) + ...
    %   tmp;
    old=final(region{:});
    mask=tmp>0.05;
    old(mask) = old(mask) * transparency;
    final(region{:}) = old + tmp;
  end
  allsubj=[final ; allsubj]; %  because 'down is up' on the images axis!
  if 0 % draw individually
    final=final/max(final(:));
    imagesc(final)
    set(gca,'ydir','normal')
    axis off
    drawnow
  end
end
if 1 % global image
  if WHOLE_FIGURE
    clf
  end
  if AMPLIFY_LESIONS
    allsubj=allsubj/max(allsubj(:));
  else
    allsubj=max(0,min(1,allsubj));
  end
  if INVERT % INVERT ?
    allsubj=1-allsubj;
  end
  image(allsubj);
  set(gca,'ydir','normal')
  axis off
  if WHOLE_FIGURE
    set(gca,'Position',[0 0 1 1])
  end
  if SHOW_NAMES % PUT NAMES ON?
    ht=size(final,1);
    textcolour={'w','k'};
    for i=1:length(sub)
      nchars = length(pf_names{sub(i)});
      text(800-nchars*5, ht*(length(sub)-i+0.8), pf_names(sub(i)), 'Color',textcolour{INVERT+1});
    end
  end
  if COLORBAR
    h=colorbar('east'); cm=colormap;
    % colour bar labels according to range of I
    set(h, 'yticklabel', arrayfun(@(x)sprintf('%0.3g',x), ...
      get(h,'ytick') *  ( nanmax(flat(I))-nanmin(flat(I)) )/ size(cm,1) + nanmin(flat(I)) ...
      , 'uniformoutput',0)   ); 
  end
end

return
%% EXAMPLE
% draw all files in dir?
names=arrayfun(@(x)x.name, dir('regress_t_*.nii.gz'), 'uniform',0);
files=cellfun(@(x)load_nii(x),names);
tmp = cat(4,files.img); 
if 0, tmp = tmp.*(tmp>2.4); % threshold positive t score
else  tmp =-tmp.*(tmp<-2.4); % or negative
end
lightbox_images_for_thesis(tmp, names);
