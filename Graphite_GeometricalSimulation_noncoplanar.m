% EXAMPLE SCRIPT TO GENERATE DIFFRACTION PATTERN FROM NONCOPLANAR DIFFRACTION GEOMETRY

% CLEAN WORKSPACE
home
%clear
%close all

% PATH TO DIFFWIZ LIBRARY
addpath(genpath('C:\Users\maher\Google drive\DIFFWIZ\'))
addpath(genpath('C:\Users\Cosmin\Desktop\Grand Diffraction Master\'))

% DEFINE LATTICE OR LOAD LATTICE FROM EXISTING STRUCTURE LIBRARY
addpath(genpath('C:\Users\Cosmin\Desktop\Cr2AlC'))
load Cr2AlC.mat
%load graphite.mat

% DEFINE DIRECTION OF CRYSTAL NORMAL
Lattice.Normal = [0 0 1]; 

% DEFINE X-RAYS
Probe.Type = 'xrays';
Probe.Energy = 8048.3; % [eV]
Probe.DiffractionGeometry = 'noncoplanar';
Probe.psi = 0.2; % This's the grazing angle

% DEFINE DETECTOR
Detector.Shape = 'square';
Detector.Size = 50; % diameter or length in mm
Detector.SpotFWHMx = 2; % Diffraction spot size on detector
Detector.SpotFWHMy = 2;
Detector.DistanceToSample = 50; % Sample-detector distance in mm
Detector.Offset = [0 25]; % offset from center of beam in mm

% MILLER INDICES TO LOOP OVER 
hkl = 0:10;

% MAIN FUNCTION
I = GeometricalSimulation1(Lattice, Probe, Detector, hkl, 1);







