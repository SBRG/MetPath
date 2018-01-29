function list=findRxnsFromTwoMets(model, met1, met2)

tmp1=findRxnsFromMets(model, met1);tmp2=findRxnsFromMets(model, met2); list=intersect(tmp1, tmp2);
printRxnFormula(model, list);
end
