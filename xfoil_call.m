function [xf] = xfoil_call(foil_name, Re, Mach, Ncrit, Iter, lower_AOA, upper_AOA, filter_val, propeller_name, varargin)
%xfoil_call A MATLAB function to call XFOIL and run a polar analysis.
%   This function automates the process of creating an XFOIL input file,
%   running XFOIL, and then reading the resulting polar data file. It is
%   designed to be used within a larger aerodynamics analysis script.
%
%   Inputs:
%   foil_name       - String with the name of the airfoil file (e.g., 'naca2412.dat')
%                     or the name of a NACA airfoil (e.g., 'NACA 2412').
%   Re              - Reynolds number for the analysis.
%   Mach            - Mach number for the analysis.
%   Ncrit           - Ncrit value for e^n transition criterion.
%   Iter            - Maximum number of iterations for the solver.
%   lower_AOA       - The starting angle of attack in degrees.
%   upper_AOA       - The ending angle of attack in degrees.
%   filter_val      - The filtering value to apply to the airfoil coordinates.
%   propeller_name  - Name of the propeller, used for creating a unique directory.
%   varargin        - Optional arguments. The first optional argument can be
%                     a boolean to control polar plotting.
%
%   Outputs:
%   xf              - An XFoil object containing the analysis results.

% Path to the XFOIL executable
current_dir = pwd;
cd XFOIL\
addpath(['Airfoil_Data\',propeller_name])

% Check for optional plot argument, which is the 10th input
plot_polar_check = 0; % Default: do not plot
if nargin > 9
    plot_polar_check = varargin{1};
end

%% Create a new instance of the XFOIL class, and set some properties
xf = XFOIL;
xf.KeepFiles = false; % Set it to true to keep all intermediate files created (Airfoil, Polars, ...)
xf.Visible = true;    % Set it to false to hide XFOIL plotting window

%% Create a NACA 5-series airfoil
xf.Airfoil = Airfoil(foil_name + ".dat");
% xf.Airfoil =  Airfoil.createNACA4('0012');
% To load an existing airfoil, use >> xf.Airfoil =  Airfoil('naca0012.dat');

%% Setup the action list

%Add five filtering steps to smooth the airfoil coordinates and help convergence
xf.addFiltering(filter_val);

%Switch to OPER mode, and set Reynolds = 3E7, Mach = 0.1
xf.addOperation(Re, Mach, Ncrit);

%Set maximum number of iterations
xf.addIter(Iter)

%Initializate the calculations
xf.addAlpha(lower_AOA,true);

%Create a new polar
xf.addPolarFile(foil_name + "_Polar.txt");

% %Calculate a sequence of angle of attack, from 0 to 25 degrees, step size of 0.1 degrees
xf.addAlpha(lower_AOA:0.5:upper_AOA);

%Close the polar file
xf.addClosePolarFile;

%And finally add the action to quit XFOIL
xf.addQuit;

%% Now we're ready to run XFOIL
xf.run
disp('Running XFOIL, please wait...')

%% Wait up to 100 seconds for it to finish... 
%It is possible to run more than one XFOIL instance at the same time
finished = xf.wait(100); 

% Read the polar data from the output file
try
    if finished
   disp('XFOIL analysis finished.')
   xf.readPolars;
    if plot_polar_check == 1
        xf.plotPolar(1); % Plot the polar if requested
    end
    else
    xf.kill;
    end
catch ME
    xf.kill;
    warning('Could not read polar file %s. It might be empty. Error: %s', polar_file, ME.message);
end

% Return to the original directory
cd(current_dir);

end