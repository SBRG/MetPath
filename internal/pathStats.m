function [pValLow, pValHigh] = pathStats(score, pathWeights, foldChangeRxnData, numPerm)

%Permutation test of how extreme a particular pathway score is given the
%data distribution

scores = zeros(numPerm,1);
for i = 1:numPerm
    curWeights = pathWeights;
    curInds = randperm(length(foldChangeRxnData));
    curFolds = foldChangeRxnData(curInds(1:length(find(curWeights))));
    curScoreProd = sum(curWeights.*curFolds');
    scores(i) = curScoreProd;
end

pValLow = length(find(score>=scores))/numPerm;
pValHigh = length(find(score<=scores))/numPerm;