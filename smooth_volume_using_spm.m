function I_s = smooth_volume_using_spm(I, SMOOTH)
% I_s = smooth_volume_using_spm(I, SMOOTH)
% calls spm_smooth with the volume.
% if 4D input, smooth each volume separately.
%   SMOOTH = voxels FWHM  
%   I ( x,y,z ) or I ( x,y,z,t )
% format of image is given by '01_Bassa.nii.gz' in current directory
% function is NOT thread safe as it uses temp files
% sgm 2014

% needs spm
% addpath /media/manohar/Programs/Science/MATLAB/R2012b/spm12b
if ~exist('SMOOTH', 'var'), SMOOTH=5; end

if exist('spm_smooth','file')

  nii=load_nii([getExperimentBase '/matlib/MNI152_T1_2mm.nii.gz']);
  for i=1:size(I,4)
    % spm smooths to file
    nii.img=I(:,:,:,i);
    save_nii(nii,'test.nii')
    spm_smooth('test.nii', 's_test.nii', [1 1 1]*SMOOTH)
    s_nii=load_nii('s_test.nii'); % reload file
    I_s(:,:,:,i)=s_nii.img; % I_s is the smoothed version
  end
else
  smth = @(x) convn( x, gaussKern([1 1 1]*SMOOTH)  ,'same');
  I_s = smth(I);
end
