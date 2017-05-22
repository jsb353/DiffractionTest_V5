% Relativistic Hartree-fock X-ray and electron scattering Factors, P.A.
% Doyle and P.S. Turner, Acta Cryst. (1968). A24, 390

function f = ScatteringFactor(s,Z,varargin)

ind = 1; % electrons

if nargin<2
    error('Insufficient number of inputs. You should supply the Scattering vector and the Atomic number. e.g. ScatteringFactor(1.02, 6)');
end

if size(varargin,2)==1
    if strcmp(varargin{1},'xrays')==1
        ind = 2;
    end
end

s = s/2; % Doyleand Turner define s as sin(theta)/lambda

m  = 9.109389754e-31; % electron mass
c = 299792458; % speed of light
q = 1.60217662e-19; % elementary  charge

 % a1 b1 a2 b2 a3 b3 a4 b4 c
Params{1}= zeros(92,9); % electrons
Params{2}= zeros(92,9); % x-rays

% Carbon
Params{2}(6,:)  = [2.3100 20.8439 1.0200 10.2075 1.5886 0.5687 0.8650 51.6512 0.2156]; 

% Oxygen
Params{2}(8,:)  = [3.0485 13.2771 2.2868 5.7011 1.5463 0.3239 0.8670 32.9089 0.2508]; 

% Silicon
Params{1}(14,:) = [2.1293 57.7748 2.5333 16.4756 0.8349 2.8796 0.3216 0.3860 0]; % a1 b1 a2 b2 a3 b3 a4 b4 c

% Vanadium
Params{2}(23,:)  = [10.297 6.866 7.351 0.438 2.070 26.894 2.057 102.478 1.220]; 

% Germanium
Params{2}(32,:)  = [16.0816 2.8509 6.3747 0.2516 3.7068 11.4468 3.6830 54.7625 2.1313]; % needs rechecking from paper copy

% Au
Params{1}(79,:) = [2.3880 42.8656 4.2259 9.7430 2.6886 2.2641 1.2551 0.3067 0]; % a1 b1 a2 b2 a3 b3 a4 b4

f = ones(size(s,1),size(Z,2)); 
for j=1:size(Z,2)
    f(:,j) = f(:,j)*Params{ind}(Z(j),9); % c
    for i=1:2:7
        f(:,j) = f(:,j) + Params{ind}(Z(j),i)*exp(-Params{ind}(Z(j),i+1)*s.^2);
    end
end

if ind==1
   % f=f*((1000*55*q)/(m*c*c) + 1); % 55 kev
end
% See http://henke.lbl.gov/optical_constants/asf.html. Potentially useful.