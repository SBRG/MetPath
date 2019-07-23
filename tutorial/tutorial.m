%% MetPath Tutorial Start

% Not necessary to run clear all but can help clarity to begin 
% with an empty Workspace
clear

%% Set up paths

% MetPath uses functions from COBRA Toolbox package
% (https://opencobra.github.io/cobratoolbox/stable/).
% Its installation is mandatory.
% If you have not yet run initCobraToolbox, do so here (otherwise skip)
% run initCobraToolbox

%% Set current directory

% Change the current directory to the MetPath folder
cd('C:\MetPath')

%% Load tutorial variables


% In order to use MetPath toolbox we need:
% - a struct object (named exprData) with 
% (1) exprData.aerobic - expression data (double) of condition 1,
% (2) exprData.anaerobic - expression data (double) of condition 2,
% (3) exprData.genes - gene names.
% modelAna - iJO1366 with exchange reactions set to anaerobic growth on glc
% modelStd - iJO1366 with exchange reactions set to aerobic growth on glc
load('data\tutorialData\tutorialStart');

% We also need currency and cofactor pairs to use in the tiered pathway
% definition.
load('data\curCofPairs');

%If pre-determined cofactor sets are not desired, run this set of lines:
% cofactorPairs = {};
% compartments = {};
% currencyPairs = {};


%% Set up solvers

% Any solver handled by the COBRA toolbox can be used
changeCobraSolver('gurobi','LP');
changeCobraSolver('gurobi','QP');
changeCobraSolver('gurobi','MILP');
changeCobraSolver('gurobi','MIQP')


%% Case 1: Calculating pathways de novo

% For the anaerobic condition, convert the model to be handled by the toolbox, 
% define the flux state and generate modelMets struct object necessary for 
% the next steps with:

biomassInd = find(strncmp('BIOMASS',modelStd.rxns,length('BIOMASS')));
%Alternatively, the user can specify a biomass index if the name does not match
%biomassInd = [];

%Can ignore inorganic metabolite handling, use a custom list, or use the
%standard list
% inorganicMets = 1; %Use default list
% inorganicMets = 0; %Ignore
% inorganicMets = {}; %Custom list (user must ensure met names are valid)
inorganicMets = {'o2', 'so2','so3','so4','nh4','no2','no3','fe2','fe3',...
'h2o','co2','co','h2o2','o2s','h2s','etha', 'no','fe3hox','3fe4s',...
'4fe4s','2fe2s', 'etoh','mobd','cu','cu2'}; %Example of a custom list

% Converting the model formatting for string matching
modelAnaAdj = convertModel(modelAna); 

%This takes in currency and cofactor pairs and creates a structure out of
%them (along with inorganic metabolites too)
[metsCurCofInorg] = setupMetClasses(currencyPairs, cofactorPairs, compartments, inorganicMets);

%This returns an active flux state of the model
%This can be substituted for a user-defined representative flux state
tolFlux = 10^-6; %The lowest non-zero flux that will not be truncated
fluxes = calculateFluxState(modelAnaAdj, tolFlux);

%This returns the network that is active
[modelAnaAdjNoBM, modelMetsAna, nonCarbonMets, fluxesRed] = getActiveNetwork(modelAnaAdj,...
    biomassInd, fluxes, inorganicMets, compartments);

%Organize Gene-Protein-Reaction (GPRs) for data mapping purposes
[parsedGPR,corrRxn] = extractGPRs(modelAnaAdjNoBM);

%Map expression data to model
fMapAna = mapGenes(modelAnaAdjNoBM, parsedGPR,corrRxn, exprData.genes, ... 
    exprData.anaerobic, exprData.aerobic);


% extract the pathways and calculate the production scores, degradation
% scores and the aggregate perturbation scores
cutoffDistance = 1;
cutoffFraction = 0.05;

%Defining different S matrices for normal, cofactor, and currency
%metabolites
sPruned = pruneMatrices(modelAnaAdjNoBM, modelMetsAna, metsCurCofInorg);

%Pathway calculation
pathsAna = metPath(modelAnaAdjNoBM, modelMetsAna, metsCurCofInorg, sPruned, cutoffDistance,cutoffFraction);

%Scoring the expression for each pathway and returning a permutation
%p-value
%cRes =  vector used to generate resultsTab, useful for other functions
numPerms = 100;
cResAna = calcRes(modelAnaAdjNoBM, modelMetsAna, fMapAna, pathsAna, numPerms);

%Collecting the results
%   resultsTab = a cell array that could be sorted
resultsTab = createResultsTab(modelMetsAna, cResAna);


% NOTE: the cutoff distance is set to 1 in order to run the tutorial
% faster. For each metabolite extract the reactions involved in the
% production and in its degradation, estimate their weightings and their
% levels. The levels are meant as the distance from the reaction directly
% involved in the production of the metabolite.

% the user can set the distance of the rxns to take in consideration during
% the pathways extraction (cutoffDistance). It is not suggested to use a
% distance too high since it leads to an increment of the time needed and
% to a loss of statistical power. A distance of 2 or 3 is a good tradeoff.

% Also since some reactions barely participate in the pathway, due to low
% flux contribution fraction, it is possible to set a threshold to filter 
% out those reactions using the cutoffFraction option


% to obtain the aggregate perturbation score (APS) in a sorted cell array we can
% use the following function:

aggregatePerturbationScoresAna = calcAggregateScores(modelMetsAna, cResAna);


%% Case 2:

%MAKE THIS SECTION INTERACT BETTER - IT'S NOT ACTUALLY A SECOND SET OF
%ANALYSES, IT'S A SECOND CONDITION WHICH WE THEN COMPARE


% In order to retrieve the subSystems perturbation score it is necessary
% compare the scores of the two differents conditions. Thus, what was done
% for the anaerobic growing conditions has to be done also for the aerobic
% growing conditions:

modelStdAdj = convertModel(modelStd);
tolFlux = 10^-6; 
fluxesStd = calculateFluxState(modelStdAdj, tolFlux);
[modelStdAdjNoBM, modelMetsStd, nonCarbonMets, fluxesRed] = getActiveNetwork(modelStdAdj,...
    biomassInd, fluxesStd, inorganicMets, compartments);
[parsedGPR,corrRxn] = extractGPRs(modelStdAdjNoBM);
fMapStd = mapGenes(modelStdAdjNoBM, parsedGPR,corrRxn, exprData.genes, ... 
    exprData.anaerobic, exprData.aerobic);
cutoffDistance = 1;
cutoffFraction = 0.00;

pathsStd = metPath(modelStdAdjNoBM, modelMetsStd, metsCurCofInorg, cutoffDistance,cutoffFraction);
cResStd = calcRes(modelStdAdjNoBM, modelMetsStd, fMapStd, pathsStd, numPerms);


%% Compare the two states %CHANGE THIS NAME
% Then we can score the subSystems Perturbation with:
subSystemsPerturbation = subSystemsScores(modelAnaAdjNoBM, cResAna, modelMetsAna,modelStd, cResStd, modelMetsStd);

% this function will calculate the overall perturbation of the subSystems in the model

% to generate a table with a comparisons in terms of perturbation score and
% used rxns we can use the comparePaths function:

[commonPaths, diffPaths] = comparePaths(modelAnaAdjNoBM,modelStdAdjNoBM,...
    cResAna,cResStd, modelMetsAna,modelMetsStd);

% % NOTE: in this case, since we extracted paths using a distance = 1 they
% will result exactly the same except for the perturbation score.


% findGenesFromPaths:  it retrieves genes involved in a pathway
% (in production reactions, degradation reactions and in th whole path)


[Pgenes, Dgenes, PDgenes] = findGenesFromPaths(cResAna, modelAnaAdjNoBM, modelMetsAna);



