close all; clear; clc;
addpath(genpath('../function/'));

%% ------------------------ Parameters --------------------------------------------
isGPU = 1;

useGPU(isGPU);

data_dir = './output/mat/';
out_name = 'LS';   % LS, KL


solve_lst = dir([data_dir, out_name, '*.mat']);  
run('PlotCmp.m');

