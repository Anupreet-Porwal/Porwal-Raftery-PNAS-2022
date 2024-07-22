The last 4 columns are the wide datasets in the study, the first 10 are tall datasets.


Prediction Summaries: 

R2: 18 methods X 14 datasets matrix containing R^2 averaged over 100-train test splits. Note that this version involved centering X and Y's using full data. We have modified that and so it might affect the results a little bit.   

phat: 18 methods X 14 datasets matrix containing model size averaged over 100-train test splits. It has EMVS model size missing, which should be imputed in the next iteration by implementing default EMVS (as suggested by EMVS vignetter by the authors)

PredCoverage:  11 methods X 14 datasets matrix containing empirical coverage at 95% interval averaged over 100-train test splits and the test dataset.

MeanIntScore:  11 methods X 14 datasets matrix containing mean interval score averaged over 100-train test splits and the test dataset.

CRPS:  11 methods X 14 datasets matrix containing CRPS averaged over 100-train test splits and the test dataset.



VarSel Summaries: 


dataset.true.mod.summary: Summary of the n,p of the datasets, number of variables in the generative model and its R^2 for all 14 datasets

Oracle RMSE: A 14 dim vector containing the Oracle RMSE of the generative model for each of the datasets, calculated as square root of mean of square(s.e.) of variables in the generative model. 

NormRMSE: 18 methods X 14 datasets matrix containing average RMSE (parameter estimation) averaged over 100 bootstrap datsets, normalised column wise by dividing by oracle RMSE.  

ParameterCoverage: 11 methods X 14 datasets matrix containing average empirical coverage at 95% level averaged over parameters and 100 bootstrapped datasets. Notice that the values are very large for wide datasets (last 4 columns) and I need to do coverage seperately for zero and non zero variables. 

ParameterIntScore: 11 methods X 14 datasets matrix containing mean interval score averaged over parameters and 100 bootstrapped datasets. 

AUPRC: 17 methods X 14 datasets matrix containing Area under precision recall curve (using the variable in generative model as true variables) averaged over 100 bootstrapped datasets. Note that there are only 17 techniques because LASSO has just one entry (since we find PIP proxy by varying over cross validation parameter so there are no two entries for Lasso - lambdamin and Lasso- Lambda.1se)

AUROC: 17 methods X 14 datasets matrix containing Area under ROC (using the variable in generative model as true variables) averaged over 100 bootstrapped datasets. Note that there are only 17 techniques because LASSO has just one entry (since we find PIP proxy by varying over cross validation parameter so there are no two entries for Lasso - lambdamin and Lasso- Lambda.1se)

comptime: 18 methods X 14 datasets matrix containing average computation time averaged over 100 bootstrap datsets.

