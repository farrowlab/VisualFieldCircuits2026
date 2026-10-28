function lor = log_odds_ratio(pU, pL)
%LOG_ODDS_RATIO Log odds ratio of two proportions, log(odds(pU) / odds(pL)).

odds_U = pU / (1 - pU);
odds_L = pL / (1 - pL);
lor = log(odds_U / odds_L);
end
