% Geometric simulation of x-ray diffraction on a surface.

home;
close all;

addpath(genpath('C:\Users\Cosmin\Desktop\Grand Diffraction Master'))
addpath(genpath('C:\Users\Cosmin\Desktop\Cr2AlC'))
% addpath(genpath('C:\Users\Cosmin\Desktop\Diffraction-master\StructureLibrary'))
% addpath(genpath('C:\Users\Cosmin\Desktop\Diffraction-master\TestScripts'))

% Sample
load NaCl.mat

% Probe    
Probe.Type = 'x-ray';
Probe.Energy = 12000; % [eV] %10000
Probe.DiffractionGeometry = 'noncoplanar';
Probe.psi = 0.4; % This is the grazing angle

% Detector
Detector.Shape = 'square'; %circle
Detector.Size = 102; % diameter or length in mm %15
Detector.SpotFWHMx = 3; % Diffraction spot size on detector % .3
Detector.SpotFWHMy = 3; %.3
Detector.DistanceToSample = 50; % Sample-detector distance in mm %200
Detector.Offset = [0 40]; % offset from center of beam in mm [0 0]

% Lattice normal of the surface
Lattice.Normal = [0 0 1];

% Miller indeces to loop over
hkl=[0:9]; %0:10

Table=GeometricSimulationNonCop(Lattice,Probe, hkl);
