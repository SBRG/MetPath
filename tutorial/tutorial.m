%% MetPath Tutorial Start

% Not necessary to run clear all but can help clarity to begin 
% with an empty Workspace
clear

%% Set up paths

% MetPath uses functions from COBRA Toolbox package
% (https://opencobra.github.io/cobratoolbox/stable/).
% Its installation is mandatory.
% If you have not yet run initCobraToolbox, do so here (otherwise skip)
run initCobraToolbox

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
%THIS WOULD DELETE THE PAIRS - MAKE SURE THIS IS HANDLED PROPERLY
cofactorPairs = {};
compartments = {};
currencyPairs = {};


%% Set up solvers

% Any solver handled by the COBRA toolbox can be used
changeCobraSolver('gurobi6', 'MIQP')
changeCobraSolver('gurobi6','LP');
changeCobraSolver('gurobi6','QP');
changeCobraSolver('gurobi6','MILP');


%%

% for anaerobic condition convert the model to be handled by the toolbox, 
% define the flux state and generate modelMets struct object necessary for 
% the next steps with:

%THIS IS STILL BROKEN WHEN IT ACTUALLY FINDS BIOMASS ENTRIES
biomassInd = find(strncmp('BIOMASS',modelStd.rxns,length('BIOMASS')));
%biomassInd = [];

allowLoops = 1;
%THIS WHOLE TUTORIAL COMBINES SNAKE CASE AND CAMEL CASE. PICK ONE BASED ON
%HOW THE FUNCTIONS ARE CALLED
%THIS SHOULD TAKE EMPTY BRACE INSTEAD OF 0
%NEED TO OFFER AN EXPLANATION OF WHAT PROVIDING AN INORGANIC METS LIST DOES
inorganicMets = 0;

% inorganicMets = {'o2', 'so2','so3','so4','nh4','no2','no3','fe2','fe3',...
% 'h2o','co2','co','h2o2','o2s','h2s','etha', 'no','fe3hox','3fe4s',...
% '4fe4s','2fe2s', 'etoh','mobd','cu','cu2'};

%ACTUALLY THE WAY WE SHOULD HANDLE BIOMASS WOULD BE TO SPLIT IT UP.
%CALCULATE THE FLUX STATE WITH BIOMASS INCLUDED, THEN SPLIT THE BIOMASS
%REACTION UP FOR PATHWAY CALCULATION. BY REMOVING THE BIOMASS, WE'RE
%IGNORING A MAJOR NON-GENE ASSOCIATED DRAIN ON MANY METABOLITES

%THIS FUNCTION IS CALLED DEFINE METS BUT IT REALLY CALCULATES FLUX AND
%CORRECTS MODEL NAMES AND A BUNCH OF UNRELATED THINGS - NEED TO RENAME OR
%SPLIT IT UP
%THIS FUNCTION CURRENTLY REMOVES BIOMASS FROM THE ORIGINAL MODEL VARIABLE
%WITHOUT RENAMING IT. AWFUL BEHAVIOR

%MAKE SURE EACH FUNCTION HAS A DESCRIPTION ABOVE IT FOR THE TUTORIAL
%(STANDARD COMMENT)
%REFACTORED FUNCTIONS START HERE
modelAnaAdj = convertModel(modelAna); %THIS FUNCTION DOESN'T SEEM NECESSARY, IT CONVERTS NAMES BUT WHY DO WE NEED THAT

%This takes in currency and cofactor pairs and creates a structure out of
%them (along with inorganic metabolites too)
[metsCurCofInorg] = setupMetClasses(currencyPairs, cofactorPairs, compartments, inorganicMets);

%This returns an active flux state of the model
mode = 3;
solFinalVals = calculateFluxState(modelAnaAdj, mode, allowLoops);

%This returns the network that is active
[modelAnaAdjNoBM, modelMetsAna, nonCarbonMets, fluxesRed] = getActiveNetwork(modelAnaAdj,...
    biomassInd, solFinalVals, inorganicMets);

%MAKE SURE THESE ARE RENAMED 
%CHECKMETS IS REMOVED TOO
% [modelAnaAdjNewNew,modelMets,otherMets] = checkMets(modelAnaAdjNew,biomassInd, inorganicMets, currencyPairs,...
%     cofactorPairs, compartments);
%REMOVE THE DEFINE METS FUNCTION
%RIGHT NOW THE FLUX STATE HAS FEWER DIMENSIONS THAN THE MODEL ITSELF
%[modelAnaMod,modelMetsAna,nonCarbonMets] = defineMets(modelAna, biomassInd, 2, ...
%    currencyPairs, cofactorPairs, compartments,allowLoops,inorganicMets);

%FIX THIS EXPLANATION TO BE MORE CLEAR
% % NOTE: biomass index must be provided, in fact, it involves too many
% % metabolites and otherwise it will be included in several pathways 
% % altering results.you can find biomass_ind by using:
% % biomassInd = strmatch('BIOMASS', modelAna.rxns);

% map the fold change of genes expression onto active reactions.
%NOT SURE WHAT THIS NOTE IS TALKING ABOUT
% % NOTE: the first expression data must refer to the model and
% % condition used for the function

[parsedGPR,corrRxn] = extractGPRs(modelAnaAdjNoBM); %IS THIS A COBRA TOOLBOX FUNCTION OR DO WE NEED TO INCLUDE IT?

%WHAT DOES FMAP STAND FOR?? TRY TO COME UP WITH A MORE CLEAR NAME
fMapAna = mapGenes(modelAnaAdjNoBM, parsedGPR,corrRxn, exprData.genes, ... 
    exprData.anaerobic, exprData.aerobic);

%RENAME metsCurCofInorg to something less specific and easier to read

%THE MODEL HAS NOT BEEN MODIFIED AND THUS DOESN"T NEED TO BE RETURNED OR
%RENAMED BY THESE FUNCTIONS. INSTEAD, NEED TO DIRECTLY RETURN RESULTS LIKE
%curCof AND PASS THOSE TO THE NEW FUNCTIONS

% extract the pathways and calculate the production scores, degradation
% scores and the aggregate perturbation scores
%THIS IS THE FUNCTION THAT HAS ANOTHER OUTPUT FORMAT THAT BREAKS LATER CODE
%TRY TO REMOVE THIS BEHAVIOR IF POSSIBLE
%WHAT DOES CRES STAND FOR? THESE NAMES DON"T MAKE SENSE
%SHOULD SEPARATE THE PATHWAY CALCULATION FROM THE EXPRESSION DATA SCORING
%I.E. FMAPANA SHOULD NOT BE PASSED TO IT, IT SHOULD BE IN A SEPARATE
%FUNCTION
%NEED TO USE GUROBI in solve_milp for the SOLVER OPTION - THIS ISN'T GREAT
%BEHAVIOR
cutoffDistance = 1;
cutoffFraction = 0.00;
paths = metPath(modelAnaAdjNoBM, modelMetsAna, metsCurCofInorg, cutoffDistance,cutoffFraction);

%Scoring the expression for each pathway and returning a permutation
%p-value
%THIS RETURNS A LOT MORE THAN JUST THE P-VALUES. SEEMS LIKE ANOTHER
%OPPORTUNITY TO SPLIT UP THE FUNCTION, SINCE THE BASIC STATS WOULD BE RUN
%ONCE BUT MAYBE YOU WANT TO TRY THE P-VALUE AT DIFFERENT LEVELS
%MAYBE MOVE THE REST OF IT TO createResultsTab?
numPerms = 1000;
cRes = calcRes(modelAnaAdjNoBM, modelMetsAna, fMapAna, paths, numPerms);

%Collecting the results
resultsTab = createResultsTab(modelMetsAna, cRes);

%WEIRD RESULTS FROM THE TUTORIAL
%succinate doesnt have a production pathway??? could that be
%currency/cofactor issues? and prpp doesnt have a consumption pathway???
%There are also metabolites where there is no 0 reaction in the pathway,
%i.e. no reaction directly consumes/produces it. This seems like it has to
%be an issue with the filtering

%CAN I SPEED IT UP? A DISTANCE OF 2 SHOULDN'T TAKE THAT LONG
% % NOTE: the cutoff distance is set to 1 in order to run the tutorial
% % faster. For each metabolite extract the reactions involved in the
% % production and in its degradation, estimate their weightings and their
% % levels. The levels are meant as the distance from the reaction directly
% % involved in the production of the metabolite.

% % the user can set the distance of the rxns to take in consideration during
% % the pathways extraction (cutoffDistance). It is not suggested to use a
% % distance too high since it leads to an increment of the time needed and
% % to a loss of statistical power. A distance of 2 or 3 is a good tradeoff.
% % Since some reactions, specially at longer distances, partecipate barely
% % in the pathway is possible to set a threshold to filter out that reactions
% % (cutoffFraction)


% to obtain the aggregate perturbation score (APS) in a sorted cell array we can
% use the following function:

aggregatePerturbationScores = APS(modelMetsAna, cResAna);


%% Second set of analyses

% In order to retrieve the subSystems perturbation score it is necessary
% compare the scores of the two differents conditions. Thus, what was done
% for the anaerobic growing conditions has to be done also for the aerobic
% growing conditions:

%DEFINE METS SHOULDN"T PRINT ANYTHING EITHER

%WHY IS MAP GENES SO SLOW??

[modelStdMod,modelMetsStd] = defineMets(modelStd, biomassInd, 2, currencyPairs, cofactorPairs, compartments,1,0);
fMapStd = mapGenes(modelStdMod, fc.genes, fc.aerobic, fc.anaerobic);
[resultsTab, cResStd] = metPath(modelStdMod, modelMets_std, fMapStd, 1, 0.05);


% Then we can score the subSystems Perturbation by

subSystemsPerturbation = subSystemsScores(modelAnaMod, cResAna, modelMetsAna,modelStd, cResStd, modelMetsStd);

% this function will predict the overall perturbation of the subSystems in the model


% generate an output compatible with escher: suppose we want study the
% cytosolic pyruvate pathway in the two condition, we can obtain an output
% (to cut and paste on a txt file and load on (escher.github.io) which
% describe the rxn scores fold change for each rxns in that pathway

% In order to use escherPaths function we need to set which is the first
% condition and the second condition:

fc.expression1 = fc.anaerobic;
fc.expression2 = fc.aerobic;

% In this case the lofFC values will be obtained by comparing anaerobic
% conditions vs aerobic conditions

% then we can run the option

%DOES THIS GENERATE A FILE THEN? WHY IS THE OUTPUT CALLED FILE? IS IT
%CALLED ANYWHERE ELSE IN THE TUTORIAL? THE OUTPUT IS AN OBNOXIOUS TEXT
%BLOCK
%THIS FUNCTION TAKES exprData IN WHICH HAS ANAEROBIC AND AEROBIC DEFINED
%SPECIFICALLY - WILL THIS WORK FOR NEW CONDITIONS???
file = escherPaths('pyr_c', 'b',modelAnaMod,exprData, cResAna,modelMetsAna,cResStd,modelMetsStd);

% alternatively we can study the pathway in a single condition by:

file = escherPaths('pyr_c', 'b',modelAnaMod,exprData, cResAna,modelMetsAna);

% % NOTE: the 'b' stands for 'both' (production and degradation), if we want study just the production or
% % degradation for a specific metabolite we can use respectvively 'p' or
% 'd'


% to generate a table with a comparisons in terms of perturbation score and
% used rxns we can use the comparePaths function:

[commonPaths, diffPaths] = comparePaths(modelAnaMod,modelStd,cResAna,cResStd, modelMetsAna,modelMetsStd);

% % NOTE: in this case, since we extracted paths using a distance = 1 they
% will result exactly the same except for the perturbation score.


% findGenesFromPaths:  it retrieves genes involved in a pathway
% (in production reactions, degradation reactions and in th whole path)

%THESE LISTS HAVE DUPLICATES
%ALSO WHAT ARE THESE OUTPUTS? THE NAMES AREN'T CLEAR - PRODUCTION AND
%DEGRADATION? I THOUGHT WE WERE CALLING THEM PRODUCTION AND CONSUMPTION?
%CHANGE D TO C AND CHECK IN THE FUNCTION TOO FOR THE SAME NAME
[Pgenes, Dgenes, PDgenes] = findGenesFromPaths(cResAna, modelAnaMod, modelMetsAna);





%% Using the universal pathway database for E. coli

load tutorialUniversalDb

% first of all the index of the genes of expression data must match those from
% universal database. It is possible to use data from chip array or from
% RNA-seq. In the latter case, if data for standard condition are not
% provided the function will use microarray data from standard conditions.
% Thus, data values will be converted in rank score values in order to
% allow a comparison of array data and rna seq data.
% the data used in this tutorial are the same used to obtain results
% described in the paper and they come from rna seq. In order to match the
% data and to rank them we can run the ollowing function:

[data_matched] = matchExpressionData(data, 1);

% % NOTE: data must contains two fields:
% % one named genes and one named vals RNAseq is a flag value
% % specifying if used data are from RNAseq experiment, standardExpression
% % can be empty (in this case will be used data from package) or user can
% % provide his own expression data.

listCond = {'acetateaerobicNH4.mat'
    'acetateanerobicNH4.mat'
    'acetateanerobicNO3.mat'
    'ala_ana.mat'
    'ala_o2.mat'
    'arg_ana.mat'
    'arg_o2.mat'
    'asn_ana.mat'
    'asn_o2.mat'
    'asp_ana.mat'
    'asp_o2.mat'
    'cys_ana.mat'
    'cys_o2.mat'
    'fumarateaerobicNH4.mat'
    'fumarateanerobicNH4.mat'
    'fumarateanerobicNO3.mat'
    'galaerobicNH4.mat'
    'galanaerobicNH4.mat'
    'galanaerobicNO3.mat'
    'glcAerobicNH4.mat'
    'glcAnaerobicNH4.mat'
    'glcAnaerobicNH4_2.mat'
    'glcAnaerobicNo3.mat'
    'gln_ana.mat'
    'gln_o2.mat'
    'glu_ana.mat'
    'glu_o2.mat'
    'gly_ana.mat'
    'gly_o2.mat'
    'glycaerobicNH4.mat'
    'glycanaerobicNH4.mat'
    'glycanaerobicNO3.mat'
    'glycolateaerobicNH4.mat'
    'glycolateanerobicNH4.mat'
    'glycolateanerobicNO3.mat'
    'his_ana.mat'
    'his_o2.mat'
    'iso_ana.mat'
    'iso_o2.mat'
    'lacAerobicNH4.mat'
    'lacAnaerobicNH4.mat'
    'lacAnaerobicNo3.mat'
    'leu_ana.mat'
    'leu_o2.mat'
    'lys_ana.mat'
    'lys_o2.mat'
    'mannoseaerobicNH4.mat'
    'mannoseanaerobicNH4.mat'
    'mannoseanaerobicNO3.mat'
    'met_ana.mat'
    'met_o2.mat'
    'phe_ana.mat'
    'phe_o2.mat'
    'pro_ana.mat'
    'pro_o2.mat'
    'succinateaerobicNH4.mat'
    'succinateanerobicNH4.mat'
    'succinateanerobicNO3.mat'
    'thr_ana.mat'
    'thr_o2.mat'
    'trp_ana.mat'
    'trp_o2.mat'
    'tyr_ana.mat'
    'tyr_o2.mat'
    'val_ana.mat'
    'val_o2.mat'};

%  after data preparation, we can use universalDb function by means of
%  ranked genes (data_matched.rank or .standardRank)
%
[pathways, perturbationScores, universal_ssScore] = universalDb(model,data_matched.genes,data_matched.rank,data_matched.standardRank, 1, listCond);

% % NOTE: exprs1 = genes values or ranked genes list from user's data
% % expression, exprs2 = genes values or ranked genes list from standard
% % conditions. ssScore = flag to score also the subSystems score, optional
% % since it drastically increase the computational time)
%
% % Usage example:
% % [pathways, perturbationScores, universal_ssScore] = universalDb(model,data_matched.genes,data_matched.vals,data_matched.standardVals, 1)

%Ok there's no reason universalDB needs to take a year. The condition files
%already have the pathways. Need to just extract those pathways and map the
%data directly

%WHERE IS THE FUNCTION THAT LOOKS AT SIMILARITY OF PATHWAYS WHEN
%GENERATING A UNIVERSAL DB TO COLLAPSE SIMILAR PATHWAYS?



