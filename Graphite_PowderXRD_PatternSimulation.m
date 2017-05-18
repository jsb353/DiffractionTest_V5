% EXAMPLE SCRIPT TO GENERATE POWDER DIFFRACTION PATTERN

% CLEAN WORKSPACE
home
clear
close all

% PATH TO LIBRARY


% DEFINE LATTICE OR LOAD LATTICE FROM EXISTING STRUCTURE LIBRARY
load graphite.mat;

% DEFINE X-RAYS
Probe.Type='x-ray';
Probe.Energy=8048.3; % [eV]

% DEFINE PARAMETERS THAT CAN BE USED IN SIMULATION
% FigNum=3;
% hkl=7;
% Threshold=1e-2;
% Separation=0.05;


Table=Generate_Intensity_2theta(Lattice, Probe);
% Generate_Intensity_2theta(Lattice,Probe, FigNum, hkl, Threshold, Separation);
