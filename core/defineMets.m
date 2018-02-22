function [model,modelMets,otherMets] = defineMets(model, biomassInd, currencyPairs, cofactorPairs, compartments, allowLoops, inorganicMets)

%THIS FUNCTION PROBABLY SHOULDN'T EVEN EXIST. IT'S A RANDOM SUBSET OF THE
%WHOLE WORKFLOW. JUST SPLIT IT INTO ITS SEPARATE FUNCTIONS AND CALL THEM
%SEPARATELY

% this function convert the model in order to be handled by the toolbox,
% define the flux state and generate modelMets struct object necessary for
% the next steps.

%WHY IS A FLUX STATE CALCULATED AS PART OF THIS FUNCTION? CONVERT MODEL
%SHOULD BE ITS OWN COMMAND IF ANYTHING. THIS IS A WEIRD COMBINATION OF
%FUNCTIONS THAT SHOULDN'T EXIST

% FBAmode permits to set how to solve the LP problem: 1-pFBA, 2-low tolerance,
% 3-minimize fluxes fixing the lb of obj, 4- minimize fluxes

% allowLoops: 1 = allow or 0 = not allow

% inorganicMets: define the list of inorganic metabolites to take in
% consideration for pathways extraction. If not present the default list
% will be used

%THIS LINE AFFECTS THE METABOLITE NAMES. THIS NEEDS TO BE AN OPTION. ALSO,
%IT REMOVES __ WHICH I DON'T THINK IT NEEDS TO DO FOR ANY REASON. ALSO A
%LOT OF OTHER THINGS ARE CHANGED THAT I DON'T THINK IT NEEDS TO BE. 
%THE FUNCTION DOES HELP TO DEAL WITH NON-STANDARD MODELS THOUGH SO IT
%SHOULD REMAIN AN OPTION
%IF I WANT TO CHECK FOR INPUTS, DO IT LIKE THIS
% ~exist('mode', 'var') | ~ismember(mode, [1 2 3 4])

%THIS FUNCTION RETURNS A MODEL WITH A FLUX VARIABLE ATTACHED. THAT SEEMS OK
%SINCE IT HELPS KEEP TRACK OF IT, AS LONG AS MULTIPLE FLUX STATES WON'T
%NEED TO BE ASSOCIATED WITH THE SAME MODEL
%THIS IS THE PART OF THE CODE THAT ALLOWS THE USER INPUT ON THE FLUX MODE
%AT THE COMMAND PROMPT. JUST REMOVE THAT ENTIRELY





