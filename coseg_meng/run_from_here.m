
clear all;
close all;

%% parameter setting
img_path='.\001\';

para_energy=zeros(1,2);
para_energy(1)=3;% the interior parameter.the larger the para_energy(1), the more power of expanding
para_energy(2)=6;% the larger the para_energy(2), the more power of shrinking
% -----------------the parameters para_energy(2) can be changed in [2,7]
ratio=[0.1,0.5,1];%for hierarchical co-segmentation. ratio = [0.3,0.5,0.7,0.9,1]...for the method in SMCB, ratio=[1].
iter=[300,300,30];%iteration number for each layer, related to ratio
is_show=1; % 0: not show, 1: show the results

%% run the co-segmentation
out_mask=run_coseg(img_path, para_energy, is_show,ratio,iter);




