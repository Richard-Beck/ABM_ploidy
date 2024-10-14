cd('/Users/4477116/Documents/projects/Ploidy_abm/optimization/')
path2root= '/Users/4477116/Documents/projects/Ploidy_abm/input/';
addpath('/Users/4477116/Documents/projects/Ploidy_abm/input/')
addpath('/Users/4477116/Documents/projects/Ploidy_abm')


% Load initial parameters from file
initial_params = readParameters('parameters.txt');

 


% Load target population from file
target_population = load('/Users/4477116/Documents/projects/Ploidy_abm/optimization/targetPop.txt');

% Define optimization parameters
alpha = 1;  % Weight for population size deviation
beta = 1;   % Weight for location deviation




% Define the optimization problem
objFunction = @(params) objectiveFun(params, target_population,alpha,beta,initial_params);

initial_guess = [initial_params.natural_growth_rate,...
                 initial_params.natural_death_rate, ...
                 initial_params.initialTumorSize,...
                 initial_params.vssl_bdary_value];
lb = [0.043, 0.033,1000,1000];
ub = [0.083, 0.73,4000,5000];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Define options for PATTERNSEARCH %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psoptions = optimoptions('patternsearch', ...
                      'Display', ...
                      'iter', ...
                      'InitialMeshSize', 1e-1, ...
                      'MeshTolerance', 1e-3, ...
                      'OutputFcn', @saveOptimizationParams);



% Perform optimization using patternsearch
[opt_params, fval] = patternsearch(@(params)objectiveFun(params, target_population, alpha, beta,initial_params), ...
                                   initial_guess, [], [], [], [], lb, ub,psoptions);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Define GENETIC ALGORITHM options%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GAoptions = optimoptions('ga', ...
                       'Display', 'iter', ...
                       'PopulationSize', 1, ...  % Adjust as needed
                       'MaxGenerations', 50, ...  % Adjust based on problem size
                       'OutputFcn', @saveOptimizationParamsGA);

% Perform optimization using Genetic Algorithm (GA)
[opt_params, fval] = ga(@(params)objectiveFun(params, target_population, alpha, beta,initial_params), ...
                        4, ... % Number of variables (2 for now,  growth_rate and death_rate)
                        [], [], [], [], lb, ub, ... % No linear constraints, use bounds
                        [], GAoptions);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Define FMINCON options        %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fmoptions = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp','OutputFcn',  @saveOptimizationParams);

% Perform optimization
[opt_params, fval] = fmincon(@(params)objectiveFun(params, target_population,alpha,beta,initial_params), initial_guess, [], [], [], [], lb, ub, [], fmoptions);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Simulated Annealing solver      %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



global paramHistory;
paramHistory = [];




saoptions = optimoptions('simulannealbnd','OutputFcn', @saveOptimizationParams);
[x_opt, fval] = simulannealbnd(@(params)objectiveWrapper(params, target_population, alpha, beta,initial_params), initial_guess, lb, ub,saoptions);


% Display the results
disp(['Optimal birth rate: ', num2str(opt_params(1))]);
disp(['Optimal death rate: ', num2str(opt_params(2))]);
disp(['Objective function value: ', num2str(fval)]);


paramNames = {'natural_growth_rate', 'natural_death_rate'};

% Create a cell array with names and values
paramCellArray = struct('natural_growth_rate', opt_params(1), 'natural_death_rate', opt_params(2));



% Write the final optimized parameters to  file
writeParameters('/Users/4477116/Documents/projects/Ploidy_abm/input/parameters.txt', paramCellArray)

 
