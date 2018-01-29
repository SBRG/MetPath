%  solver tol = 10^12



%  glucose aerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);
model.lb(strmatch('EX_o2_e', model.rxns))=-20; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\glcAerobicNH4')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')




%  glucose anaerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);
model.lb(strmatch('EX_o2_e', model.rxns))=-2; % 1/10 of uptake in model in normal conditions 

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\glcAnaerobicNH4')

clear variables


%  glucose anaerobic NH4 -2 using other data

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);
model.lb(strmatch('EX_glc_D_e', model.rxns))=-16; % https://scihub22266oqcxt.onion.link/https://doi.org/10.1016/j.cels.2016.08.013
model.lb(strmatch('EX_o2_e', model.rxns))=-2; % 1/10 of uptake in model in normal conditions 

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\glcAnaerobicNH4_2')

clear variables






% glucose anaerobic nitrate

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
%  runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_o2_e', model.rxns))=0;
model.lb(strmatch('EX_no3_e', model.rxns))=-20;     % to define
model.lb(strmatch('EX_nh4_e', model.rxns))=0;


optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');



solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\glcAnaerobicNo3')


clear variables



% lactate aerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_lac_L_e', model.rxns))=-13;    % to difine, it cant grow with a lb = -7.8 
model.lb(strmatch('EX_o2_e', model.rxns))=-20; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/


optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\lacAerobicNH4')
clear variables



% lactate anaerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_o2_e', model.rxns))=-4;     % it cant grow with less
model.lb(strmatch('EX_lac_L_e', model.rxns))=-20;   % no data

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\lacAnaerobicNH4')
clear variables

% lactate anaerobic nitrate

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_o2_e', model.rxns))=-0;
model.lb(strmatch('EX_no3_e', model.rxns))=-20;% no data
model.lb(strmatch('EX_lac_L_e', model.rxns))=-13;% no data

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\lacAnaerobicNo3')
clear variables


% glycerol aerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_glyc_e', model.rxns))=-16.6; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4% other data show that it should be-0.11 %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2258577/
model.lb(strmatch('EX_o2_e', model.rxns))=-20; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\glycaerobicNH4')
clear variables



% glycerol anaerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_glyc_e', model.rxns))=-16.6; % no data
model.lb(strmatch('EX_o2_e', model.rxns))=-2; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\glycanaerobicNH4')
clear variables




% glycerol anaerobic nitrate

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_glyc_e', model.rxns))=-10; % no data
model.lb(strmatch('EX_o2_e', model.rxns))=0; 
model.lb(strmatch('EX_no3_e', model.rxns))=-20;% no data

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\glycanaerobicNO3')
clear variables






% galactose aerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_gal_e', model.rxns))=-9.5; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model.lb(strmatch('EX_o2_e', model.rxns))=-20; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\galaerobicNH4')
clear variables



% galactose anaerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_gal_e', model.rxns))=-9.5; % no data
model.lb(strmatch('EX_o2_e', model.rxns))=-2; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\galanaerobicNH4')
clear variables




% galactose anaerobic nitrate

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_gal_e', model.rxns))=-9.5; % no data
model.lb(strmatch('EX_o2_e', model.rxns))=0; 
model.lb(strmatch('EX_no3_e', model.rxns))=-20;% no data

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\galanaerobicNO3')
clear variables







% mannose aerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_man_e', model.rxns))=-12; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model.lb(strmatch('EX_o2_e', model.rxns))=-20; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/


optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\mannoseaerobicNH4')
clear variables



% mannose anaerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_man_e', model.rxns))=-12;; % no data
model.lb(strmatch('EX_o2_e', model.rxns))=-2; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\mannoseanaerobicNH4')
clear variables




% mannose anaerobic nitrate

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_man_e', model.rxns))=-12;; % no data
model.lb(strmatch('EX_o2_e', model.rxns))=0; 
model.lb(strmatch('EX_no3_e', model.rxns))=-20;% no data

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\mannoseanaerobicNO3')
clear variables



% acetate  aerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_ac_e', model.rxns))=-7.3000; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model.lb(strmatch('EX_o2_e', model.rxns))=-20; 
% model.lb(strmatch('EX_no3_e', model.rxns))=-20;

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\acetateaerobicNH4')
clear variables



% acetate  anerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_ac_e', model.rxns))=-7.3000; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model.lb(strmatch('EX_o2_e', model.rxns))=-5; 
% model.lb(strmatch('EX_no3_e', model.rxns))=-20;

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\acetateanerobicNH4')
clear variables


% acetate  anerobic NO3

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_ac_e', model.rxns))=-7.3000; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model.lb(strmatch('EX_o2_e', model.rxns))=-0; 
 model.lb(strmatch('EX_no3_e', model.rxns))=-20;

 optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');
 
 
solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\acetateanerobicNO3')
clear variables



% fumarate  aerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_fum_e', model.rxns))=-8.2; % almost the mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model.lb(strmatch('EX_o2_e', model.rxns))=-20; 
% model.lb(strmatch('EX_no3_e', model.rxns))=-20;

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');

solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\fumarateaerobicNH4')
clear variables



% fumarate  anerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_fum_e', model.rxns))=-7.3000; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model.lb(strmatch('EX_o2_e', model.rxns))=-4; 
% model.lb(strmatch('EX_no3_e', model.rxns))=-20;

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\fumarateanerobicNH4')
clear variables


% fumarate  anerobic NO3

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;

model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_fum_e', model.rxns))=-7.3000; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model.lb(strmatch('EX_o2_e', model.rxns))=-0; 
 model.lb(strmatch('EX_no3_e', model.rxns))=-20;

 optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');
 
solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\fumarateanerobicNO3')
clear variables


% succinate  aerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_succ_e', model.rxns))=-9.5; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model.lb(strmatch('EX_o2_e', model.rxns))=-20; 
% model.lb(strmatch('EX_no3_e', model.rxns))=-20;

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\succinateaerobicNH4')
clear variables



% succinate  anerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;

model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_succ_e', model.rxns))=-9.5; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model.lb(strmatch('EX_o2_e', model.rxns))=-6; 
% model.lb(strmatch('EX_no3_e', model.rxns))=-20;

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\succinateanerobicNH4')
clear variables


% succinate  anerobic NO3

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_succ_e', model.rxns))=-7.3000; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model.lb(strmatch('EX_o2_e', model.rxns))=-0; 
 model.lb(strmatch('EX_no3_e', model.rxns))=-20;

 
 optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');
 
solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\succinateanerobicNO3')
clear variables


% glycolate  aerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_glyclt_e', model.rxns))=-11.5; % lower level starting from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model.lb(strmatch('EX_o2_e', model.rxns))=-20;
% model.lb(strmatch('EX_no3_e', model.rxns))=-20;

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\glycolateaerobicNH4')
clear variables



% glycolate  anerobic NH4

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_glyclt_e', model.rxns))=-11.5;
model.lb(strmatch('EX_o2_e', model.rxns))=-5; 
% model.lb(strmatch('EX_no3_e', model.rxns))=-20;

optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');



solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\glycolateanerobicNH4')
clear variables


% glycolate  anerobic NO3

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
% runMetChange
model=iJO1366;


model = convertModel(model);

model.lb(strmatch('EX_glc_D_e', model.rxns))=0;
model.lb(strmatch('EX_glyclt_e', model.rxns))=-8; 
model.lb(strmatch('EX_o2_e', model.rxns))=-0; 
model.lb(strmatch('EX_no3_e', model.rxns))=-20;


optimizeCbModel(model, 'max', 10-6);
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f; % mean https://www.ncbi.nlm.nih.gov/pmc/articles/PMC294601/?page=4
model=changeObjective(model, 'ATPM');



solFinalVals = calculateFluxState(model,2,1);    
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\glycolateanerobicNO3')
clear variables
















%  ala 

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_ala_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\ala_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_ala_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\ala_ana')

clear variables


% arg

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_arg_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\arg_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_arg_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\arg_ana')

clear variables




% asn

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_asn_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\asn_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_asn_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\asn_ana')

clear variables



% asp


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_asp_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\asp_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_asp_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\asp_ana')

clear variables



% cys

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_cys_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\cys_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_cys_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\cys_ana')

clear variables


% glu

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glu_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\glu_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_glu_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\glu_ana')

clear variables



% gln

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_gln_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\gln_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_gln_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\gln_ana')

clear variables

% gly

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_gly_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\gly_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_gly_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\gly_ana')

clear variables



% his

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_his_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\his_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_his_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\his_ana')

clear variables




% iso


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_iso_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\iso_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_iso_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\iso_ana')

clear variables




% leu

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_leu_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\leu_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_leu_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\leu_ana')

clear variables





% lys

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_lys_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\lys_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_lys_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\lys_ana')

clear variables





% met

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_met_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\met_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_met_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\met_ana')

clear variables




% phe


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_phe_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\phe_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_phe_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\phe_ana')

clear variables





% pro

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_pro_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\pro_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_pro_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\pro_ana')

clear variables








% thr


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_thr_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\thr_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_thr_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\thr_ana')

clear variables








% trp

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_trp_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\trp_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_trp_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\trp_ana')

clear variables









% tyr

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_tyr_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\tyr_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_tyr_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\tyr_ana')

clear variables









% val

load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_val_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -20;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\val_o2')

clear variables


load('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\iJO1366.mat')
model = iJO1366;
model = convertModel(model);

model.lb(strmatch('EX_val_L_e', model.rxns)) = -0.5;
model.lb(strmatch('EX_o2_e', model.rxns)) = -2;

model = changeObjective(model, 'BIOMASS_Ec_iJO1366_core_53p95M');
optimizeCbModel(model, 'max', 10-6)
model.lb(strmatch('BIOMASS_Ec_iJO1366_core_53p95M', model.rxns)) = ans.f;
model = changeObjective(model, 'ATPM');


solFinalVals = calculateFluxState(model,2,1);
biomass_ind = strmatch('BIOMASS', model.rxns);
[model,modelMets,solFinalVals] = newcheckMets(solFinalVals,model,biomass_ind);
[curCof] = newcurCorPairs(modelMets);
s = pruneMatrices5(model, modelMets, curCof);
calcPaths_dist2001 = newcalcPathways(model,modelMets,curCof,s,2,0.01)
[Ppaths Dpaths]= extractPaths(model,modelMets, calcPaths_dist2001);

save('C:\Users\Palsson Lab\Dropbox (Personal)\MetChange 2.0\Gianluca\metchange2\ecoli\different condition paths extraction\val_ana')

clear variables

























