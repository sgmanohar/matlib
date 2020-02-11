function folder = getExperimentBase()
% find the main folder of Sanjay's Experiment directory

list = {
  's:/Experiment';
  '~/Documents/Experiment'; % WORK
  '~/Experiment'; % HOME
  '/media/manohar/Data/Experiment'
  '/media/DATA/smanohar/Documents/Experiment'
  'e:/smanohar/Documents/Experiment'
  'C:\Users\Sanjay\Documents\Experiment'
  'C:\Users\cogneuro\Documents\Sanjay\Experiment'
  'D:\Experiment'
  };

for i=1:length(list)
  if exist(list{i},'dir')
    break;
  end
end
folder = list{i};