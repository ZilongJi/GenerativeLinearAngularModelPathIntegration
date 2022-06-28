%% Main entry point for the Vector Addiction Modelling (VAM) project
% This script loads the data the from the Howett 2019., Brain paper. It
% prepare the variables and runs the modelling described in the current
% paper. Each figure will have a separate script that can be run after this
% script has been correctly executed. For additional dependecies please
% refer to the README file.

%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default'); %for code reproducibility

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/HowettBrain2019_Dataset.mat');
savefolder = pwd + "/Output/";