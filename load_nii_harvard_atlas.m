function [harv2a regions] = load_nii_harvard_atlas
%% ROI with  harvard-oxford atlas
fsldir = '/usr/share/fsl';
if ~exist(fsldir)
  fsldir = '/usr/local/fsl';
end
atlasses = [fsldir '/data/atlases/HarvardOxford/'];
% this will throw an error if we can't find the atlas files. 
% may need to change the above line!
harv=load_nii([atlasses 'HarvardOxford-cort-maxprob-thr25-2mm.nii.gz']);
  
harv=harv.img;
% list of regions from atlas
regions = {
'Frontal Pole'
'Insular Cortex'
'Superior Frontal Gyrus'
'Middle Frontal Gyrus'
'Inferior Frontal Gyrus, pars triangularis'
'Inferior Frontal Gyrus, pars opercularis'
'Precentral Gyrus'
'Temporal Pole'
'Superior Temporal Gyrus, anterior division'
'Superior Temporal Gyrus, posterior division'
'Middle Temporal Gyrus, anterior division'
'Middle Temporal Gyrus, posterior division'
'Middle Temporal Gyrus, temporooccipital part'
'Inferior Temporal Gyrus, anterior division'
'Inferior Temporal Gyrus, posterior division'
'Inferior Temporal Gyrus, temporooccipital part'
'Postcentral Gyrus'
'Superior Parietal Lobule'
'Supramarginal Gyrus, anterior division'
'Supramarginal Gyrus, posterior division'
'Angular Gyrus'
'Lateral Occipital Cortex, superior division'
'Lateral Occipital Cortex, inferior division'
'Intracalcarine Cortex'
'Frontal Medial Cortex'
'Juxtapositional Lobule Cortex (formerly Supplementary Motor Cortex)'
'Subcallosal Cortex'
'Paracingulate Gyrus'
'Cingulate Gyrus, anterior division'
'Cingulate Gyrus, posterior division'
'Precuneous Cortex'
'Cuneal Cortex'
'Frontal Orbital Cortex'
'Parahippocampal Gyrus, anterior division'
'Parahippocampal Gyrus, posterior division'
'Lingual Gyrus'
'Temporal Fusiform Cortex, anterior division'
'Temporal Fusiform Cortex, posterior division'
'Temporal Occipital Fusiform Cortex'
'Occipital Fusiform Gyrus'
'Frontal Operculum Cortex'
'Central Opercular Cortex'
'Parietal Operculum Cortex'
'Planum Polare'
'Heschl Gyrus'
'Planum Temporale'
'Supracalcarine Cortex'
'Occipital Pole'

};
% load subcortical atlas too
harvs=load_nii([atlasses 'HarvardOxford-sub-maxprob-thr25-2mm.nii.gz']); harvs=harvs.img;
regionss = {
'Left Cerebral White Matter'
'Left Cerebral Cortex '
'Left Lateral Ventrical'
'Left Thalamus'
'Left Caudate'
'Left Putamen'
'Left Pallidum'
'Brain-Stem'
'Left Hippocampus'
'Left Amygdala'
'Left Accumbens'
'Right Cerebral White Matter'
'Right Cerebral Cortex '
'Right Lateral Ventricle'
'Right Thalamus'
'Right Caudate'
'Right Putamen'
'Right Pallidum'
'Right Hippocampus'
'Right Amygdala'
'Right Accumbens'
};

% compile both into a single list of regions
harv = single(harv) + (harvs>0 & harv==0) .* (48 + harvs); % subcortical should start at 49
regions = [regions; regionss];

% and now load the full  probabilistic atlas
harv2=load_nii([atlasses 'HarvardOxford-cort-prob-2mm.nii.gz']); harv2=harv2.img;
harv2s=load_nii([atlasses 'HarvardOxford-sub-prob-2mm.nii.gz']); harv2s=harv2s.img;
harv2a=cat(4, harv2,harv2s); % HARV2S ( X, Y, Z, ROI )
