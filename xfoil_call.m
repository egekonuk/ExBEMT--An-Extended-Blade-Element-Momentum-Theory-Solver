function [xf] = xfoil_call(foil,Re,Mach,NCrit,N,alphaL,alphaU,filt,~,Prop_Path, plot_polars)
%Switch dir
 % addpath('XFOIL\')
 addpath(genpath('XFOIL\'));
 %cd XFOIL\
 
% addpath(['Airfoil_Data\',Sel_Prop]) % This line is no longer needed and has been removed.

%% Create a new instance of the XFOIL class, and set some properties
xf = XFOIL;
xf.KeepFiles = false; % Set it to true to keep all intermediate files created (Airfoil, Polars, ...)
xf.Visible = true;   % Set it to false to hide XFOIL plotting window

%% MODIFIED: Load an airfoil using a full path from the propeller package
% Construct the full path to the airfoil file inside the propeller's 'foils' folder

% airfoil_path = fullfile('..', 'Propeller Packages', Sel_Prop, 'foils', [foil, '.dat']);

airfoil_path = fullfile(Prop_Path, 'foils', [foil, '.dat']);

% polar_filename = fullfile('..', 'Propeller Packages', Sel_Prop, [foil, '_Polar.txt']);

polar_filename = fullfile([foil, '_Polar.txt']);
if ~exist(airfoil_path, 'file')
    error('Airfoil file not found at: %s', airfoil_path);
end
xf.Airfoil = Airfoil(airfoil_path);

%% Setup the action list

%Add five filtering steps to smooth the airfoil coordinates and help convergence
xf.addFiltering(filt);

%Switch to OPER mode, and set Reynolds, Mach, and NCrit
xf.addOperation(Re, Mach, NCrit);

%Set maximum number of iterations
xf.addIter(N)

%Initializate the calculations
xf.addAlpha(0,true);

%Create a new polar file
xf.addPolarFile(polar_filename);

% %Calculate a sequence of angle of attack, from alphaL to alphaU, step size of 1 degrees
xf.addAlpha(alphaL:1:alphaU);

%Close the polar file
xf.addClosePolarFile;

%And finally add the action to quit XFOIL
xf.addQuit;

%% Now we're ready to run XFOIL
xf.run
disp(['Running XFOIL for ', foil, ', please wait...'])

%% Wait up to 100 seconds for it to finish...
finished = xf.wait(50);

%% If successfull, read and plot the polar
if finished
    % disp('XFOIL analysis finished.') % Reduced console output
    xf.readPolars;
    % To prevent a figure popping up for every XFOIL run, the plot call is commented out.
    % figure(6)
    if plot_polars==true, figure(5); xf.plotPolar(1); else,end
    % xf.plotPolar(1);
else
    xf.kill;
    warning('XFOIL analysis for %s did not finish.', foil);
end
% cd ..