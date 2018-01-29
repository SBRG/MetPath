function [data_matched] = matchExpressionData(data, RNASeq, standardExpression)

% this simple function match the genes from experimental data to those from
% standard conditions. data must be a struct object with genes (gene
% genes - Blattner number) and expression values.
% RNASeq = flag. If 1 convert experimental data and standard expression data
% from the universal database in a ranked list in order to allow comparison
% of RNA-seq data to microarray expression data
% User can use own expression data set for standard condition or use the
% data set provided with the package.
% it returns a struct object with genes, vals = values from experimental
% data, rank = ranked list of genes from experimental data, standardVals =
% values from standard conditions, standardRank = ranked genes list from
% standard conditions.


if ~exist('RNASeq', 'var')
    RNASeq = 0;
end
if ~exist('standardExpression', 'var')
    load('standardExpression');
end

% collapsing values
tmp1 = collapseMean([data.genes, num2cell(data.vals)]);
data.genes = unique(data.genes);
data.vals = cell2mat(tmp1(:,2));

tmp2 = collapseMean([standardExpression.genes, num2cell(standardExpression.vals)]);
standardExpression.genes = unique(standardExpression.genes);
standardExpression.vals = cell2mat(tmp2(:,2));


% matching genes
lind = ismember(data.genes, standardExpression.genes);
data.genes = data.genes(lind);
data.vals = data.vals(lind);

lind = ismember(standardExpression.genes, data.genes);
standardExpression.genes = standardExpression.genes(lind);
standardExpression.vals = standardExpression.vals(lind);


if RNASeq == 1
   [data.vals, sortInd] = sort(data.vals, 'descend');
   data.genes = data.genes(sortInd);
   data.rank(:,1) = 1:numel(data.vals);
   data.rank = flipud(data.rank);
   
   [standardExpression.vals, sortInd] = sort(standardExpression.vals, 'descend');
   standardExpression.genes = standardExpression.genes(sortInd);
   standardExpression.rank(:,1) = 1:numel(standardExpression.vals);
   standardExpression.rank = flipud(standardExpression.rank);
 
   ind = match(standardExpression.genes, data.genes);
   data.genes = data.genes(ind);
   data.vals = data.vals(ind);
   data.rank = data.rank(ind);
   
else
    ind = match(standardExpression.genes, data.genes);
    data.genes = data.genes(ind);
    data.vals = data.vals(ind);
end

data_matched = data;
data_matched.standardVals = standardExpression.vals;
data_matched.standardRank = standardExpression.rank;

