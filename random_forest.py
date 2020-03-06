#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 14:04:22 2019

@author: jonathansong
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 13:38:38 2018
@author: jonathansong
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import model_selection
from sklearn.model_selection import cross_val_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)
from scipy import interp
import matplotlib.pyplot as plt

from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold

# from N_and_C_terminus import N_and_C_terminus

data = pd.read_csv("/Users/jonathansong/Research_Folder/clean_data.tsv", sep="\t")

# data = N_and_C_terminus(data_input=data)
# data = aaComp(data_input=data)

#del_vars = []
#data_vars = data.columns.values.tolist()
#data_final_vars = [i for i in data_vars if i not in del_vars]
# data_final = data[data_final_vars]

to_keep = [
 'Validation_outcome',
 'Tiny.mol.p',
 'Small.mol.p',
 'Aliphatic.mol.p',
 'Aromatic.mol.p',
 'Nonpolar.mol.p',
 'Polar.mol.p',
 'Charged.mol.p',
 'Basic.mol.p',
 'Acidic.mol.p',
 'Sulfur.mol.p',
 'N.A',
 'membrane',
 'cytosol',
 'nucleus',
 'secretory.pathway',
 'mitochondrion',
 'extracellular',
 'Peptide.length',
 'Molecular.weight',
 'Net.charge',
 'Isoelectric.point',
 'Hydrophobicity',
 'Aliphatic.index',
 'Instability.index',
 'Boman.index']

data_final = data[to_keep]
data_final_vars = data_final.columns.values.tolist()

y_final_vars=["Validation_outcome"]
X_final_vars=[i for i in data_final_vars if i not in y_final_vars]

X = data_final[X_final_vars]
y = data_final[y_final_vars]

X = X.values
y = y.values


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=None)

from sklearn.ensemble import RandomForestClassifier
clf = RandomForestClassifier(n_estimators=100, class_weight='balanced')
clf.fit(X_train, y_train)

y_pred = clf.predict(X_test)
# y_score = clf.predict_proba(X_test)
weights = clf.feature_importances_
for i in range(len(weights)):
    if weights[i] > 0.08:
        print(X_final_vars[i], "--", weights[i])

print('Accuracy of random forest classifier on test set: {:.2f}'.format(clf.score(X_test, y_test)))
kfold = model_selection.KFold(n_splits=5, random_state=None)
results = cross_val_score(clf, X_train, np.ravel(y_train), cv=kfold)
print("5-fold cross validation average accuracy: %.3f" % (results.mean()))

matrix = confusion_matrix(y_test, y_pred)
print(matrix)
print(classification_report(y_test, y_pred))

#from sklearn.metrics import roc_auc_score
#y_test = y_test.astype(int)
#y = y.as_matrix()
#y_pred = y_pred.astype(int)
#logit_roc_auc = roc_auc_score(y_test, y_score)
#print("ROC AUC:", round(logit_roc_auc, 3))

#from sklearn.metrics import average_precision_score
#logit_pr_auc = average_precision_score(y_test, y_pred)
#print("PR AUC:", round(logit_pr_auc, 3))

#from sklearn.metrics import roc_curve
#import matplotlib.pyplot as plt
#fpr, tpr, thresholds = roc_curve(y, y_pred)
#fpr, tpr, thresholds = roc_curve(y_test, y_score)
#plt.figure()
#plt.plot([0, 1], [0, 1], 'r--')
#plt.axis([0, 1, 0, 1])
#plt.xlabel('False Positive Rate')
#plt.ylabel('True Positive Rate')
#plt.show()


          
cv = StratifiedKFold(n_splits=5)
   
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)

i = 0
for train, test in cv.split(X, y):
    probas_ = clf.fit(X[train], y[train]).predict_proba(X[test])
    # Compute ROC curve and area the curve
    fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
    tprs.append(interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
    plt.plot(fpr, tpr, lw=1, alpha=0.3,
             label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))

    i += 1
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
         label='Chance', alpha=.8)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
plt.plot(mean_fpr, mean_tpr, color='b',
         label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
         lw=2, alpha=.8)

std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                 label=r'$\pm$ 1 std. dev.')

plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC of Forest')
plt.legend(loc="lower right")
plt.show()