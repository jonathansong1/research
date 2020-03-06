#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 17:38:04 2019

@author: jonathansong
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 18:18:04 2019

@author: jonathansong
"""

import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
import statsmodels.api as sm
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)
from N_and_C_terminus import N_and_C_terminus
from Amino_acid_composition import aaComp



data = pd.read_csv("/Users/jonathansong/Research_Folder/out_data.tsv", sep="\t")


data = N_and_C_terminus(data_input=data)
data = aaComp(data_input=data)

del_vars = ['X',
            'Identifier',
            'Patient.ID',
            'Gene.symbol',
            'Mutant.peptide',
            'Wild.type.peptide']
data_vars = data.columns.values.tolist()
data_final_vars = [i for i in data_vars if i not in del_vars]
data_final = data[data_final_vars]
 
y=["Validation.Outcome"]
X=[i for i in data_final_vars if i not in y]

logreg = LogisticRegression(class_weight='balanced')

X=data_final[X]
y=data_final["Validation.Outcome"]
             
#X = sm.add_constant(X)
logit_model=sm.Logit(y.astype(float), X.astype(float)).fit()
print(logit_model.summary())

logreg.fit(X, y)
y_pred = logreg.predict(X)
print('Accuracy of logistic regression classifier: {:.2f}'.format(logreg.score(X, y)))

matrix = confusion_matrix(y, y_pred)
print(matrix)

print(classification_report(y, y_pred))

from sklearn.metrics import roc_auc_score
y = y.astype(int)
y = y.as_matrix()
y_pred = y_pred.astype(int)
logit_roc_auc = roc_auc_score(y, y_pred)
print("ROC AUC:", round(logit_roc_auc, 3))

print("Positive predictive value:", round((matrix[1][1] / (matrix[1][1] + matrix[0][1])), 3))

from sklearn.metrics import average_precision_score
logit_pr_auc = average_precision_score(y, y_pred)
print("PR AUC:", round(logit_pr_auc, 3))

from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt
#fpr, tpr, thresholds = roc_curve(y, y_pred)
fpr, tpr, thresholds = roc_curve(y, logreg.predict_proba(X)[:,1])
plt.figure()
plt.plot(fpr, tpr, linewidth=2, label=None)
plt.plot([0, 1], [0, 1], 'r--')
plt.axis([0, 1, 0, 1])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.show()