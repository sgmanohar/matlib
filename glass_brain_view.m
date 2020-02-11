function glass_brain_view(I, varargin)
% glass brain view
%  I = 4-dimensional volume, X,Y,Z,SUBJECT
%  Dim : [] (default) - plot all 3 dimensions in 3 subplots
%        1/2/3  plot only one dimension on current axes
%  TemplateColor : [r g b] = colour for the brain template
%                  default [0.3 0.3 0.3]  = grey
%  Color      : the color of the overlay image, default magenta
%  Background : default [0 0 0]
%  Transform  : applied to the summed values of I, before multiplying the
%               color. e.g. @(x)log(x) for log transform, or 
%                           @(x)x>0    to binarise the image
%  TemplateSlice : if you just want a single slice as the background
%                  template, specify the slice (can be integer index, or a
%                  proportion e.g. 0.5, of volume's size)
%  Boolean    : 1 means count the number of nonzero voxels at each location,
%               across volumes. 0 means just add up all the values and display.
%  you can provide a set of coordinates [x,y,z] instead of I, and I will draw
%  spheres around them.
% sgm

[ DIM col0 col1 bgcol transform tslice BOOLEAN]=parsepvpairs(...
  {'Dim', 'TemplateColor', 'Color', 'Background', 'Transform', 'TemplateSlice', 'Boolean'} , ...
  { [] , [1 1 1]*0.3, [1 0 1] , [0 0 0], @(x)x  , [] , true} , ...
    varargin{:} );
  
if size(I,2)==3 % provided coordinates
  c=I; I=zeros(91,109,91);
  for i=1:size(c,1), I(c(i,1),c(i,2),c(i,3))=1; end
  I=smooth_volume_using_spm(I,7);
  BOOLEAN=false;
end

overlay=[]; % overlay = fsub('Jansons'); selection(selection==overlay)=[];
col2=[1 0 1];     % mask colour for overlay
selection = 1:size(I,4); % select all volumes in I

% read template (where from?)
templ = load_nii([ getExperimentBase '/matlib/MNI152_T1_2mm.nii.gz']);
templ = templ.img;
if ~isrow(bgcol) bgcol=bgcol'; end;


if isempty(DIM)
  dims = [1:3];
else
  dims = DIM;
end


% composite image
for dim=dims
  if length(dims)>1
    subplot(2,2,dim);
    if dim==2, subplot(1,2,2); end
  end
  if BOOLEAN
    tmp=sq(nansum(nansum(I(:,:,:,selection)>0,dim),4))';
  else
    tmp=sq(nansum(nansum(I(:,:,:,selection),dim),4))';
  end
  tmp=transform(tmp); 
  if isempty(tslice)
    tmpl=sq(nansum(templ,dim))';
  else
    slicer = {':',':',':'};
    if tslice>1,
      slicer{dim} = tslice;
    else
      slicer{dim} = 1+floor(tslice .* size(templ, dim)); 
    end
    tmpl= double(sq(  nansum(templ(slicer{:}),dim) ))';
  end
  if any(overlay) % compose 3 items?
    tmp2=sq(nansum(nansum(I(:,:,:,overlay),dim),4))';
    tmp=bsxfun(@times, tmp,   permute(col1,[1 3 2])/max(tmp(:))  ) - ...
        bsxfun(@times, tmp2,  permute(col2,[1 3 2])/max(tmp2(:)) ) + ...
        bsxfun(@times, tmpl,  permute(col0, [1 3 2])/max(tmpl(:)) );
  else
    tmp=bsxfun(@times, tmp,   permute(col1,[1 3 2])/max(tmp(:))  ) + ...
        bsxfun(@times, tmpl,  permute(col0, [1 3 2])/max(tmpl(:)) );
  end
  tmp=bsxfun(@plus, tmp, permute(bgcol, [1 3 2] )); % add background
  tmp=max(0,min(tmp,1)); % limit to range 0-1
  imagesc(tmp); % draw
  if dim==1, set(gca, 'ydir','normal'); end
  set(gca, 'xtick',[],'ytick',[]);
  axis off
end
