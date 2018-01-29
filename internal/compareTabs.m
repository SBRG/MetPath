function comparedTab = compareTabs(agregatePerturbationScore1, aggregatePerturbationScore2)

lind_trp = ismember(agregatePerturbationScore1(:,1), aggregatePerturbationScore2(:,1));
tmp_ut_trp = agregatePerturbationScore1(lind_trp,:);
tmp_ut_trp = sortrows(tmp_ut_trp,1);



lind_std = ismember(aggregatePerturbationScore2(:,1), agregatePerturbationScore1(:,1));
tmp_ut_std = aggregatePerturbationScore2(lind_std,:);
tmp_ut_std = sortrows(tmp_ut_std,1);


ind_Trp = match(tmp_ut_trp(:,1), tmp_ut_std(:,1));

numRes = cell2mat(tmp_ut_trp(:,2:end))./cell2mat(tmp_ut_std(ind_Trp,2:end));

comparedTab = [tmp_ut_trp(:,1), num2cell(numRes)];


