import pandas as pd
import os
import sys
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix , roc_auc_score, roc_curve, precision_recall_curve, average_precision_score
from xgboost import XGBClassifier
import pickle
def pre_process(path, cut_off):
    col_names=['seqid', 'wstart','wend','ID']+list(range(80))
    df=pd.read_csv(path, sep='\t', header=None, names=col_names)
    el=df[df['ID'].str.contains('EL')]
    sl=df[df['ID'].str.contains('SL')]
    rev_idx= ['seqid','wstart','wend','ID']+list(range(79,-1,-1))
    sl_rev=sl.loc[:,rev_idx]
    sl_rev.columns=col_names
    data=pd.concat([el, sl_rev]).assign(Y= lambda x: (~x['ID'].str.contains('ref')).astype(int))
    keep=data.loc[:,list(range(80))].sum(axis=1) >= cut_off
    data_exp=data[keep]
    bal=data_exp['Y'].value_counts().min()
    data_comp=pd.concat([data_exp[data_exp['Y']==1].sample(bal, random_state=420024), data_exp[data_exp['Y']== 0].sample(bal, random_state=420024)])
    return(data_comp)


def train_model(X,Y, model=XGBClassifier(random_state=2234)):
    X_train, X_test, Y_train, Y_test = train_test_split(X,Y, test_size=.3, random_state=8976)
    model.fit(X_train, Y_train)
    Y_pred=model.predict(X_test)
    Y_prob=model.predict_proba(X_test)[:,1]
    roc_auc=roc_auc_score(Y_test, Y_prob)
    pr_auc= average_precision_score(Y_test, Y_prob)
    print('confusion matrix\n')
    print(confusion_matrix(y_pred=Y_pred,y_true=Y_test))
    print('classification_report\n')
    print(classification_report(y_pred=Y_pred,y_true=Y_test))
    print('ROC-AUC: ' + str(roc_auc))
    print('PR-AUC: ' + str(pr_auc))
    return(model)


args=sys.argv
wd=args[1]
cov_file=args[2]
trained_model_file=args[3]

os.chdir(wd)

data=pre_process(cov_file, cut_off=80)
x_cols=list(range(80))
x=data.loc[:,x_cols]
y=data['Y']
trained_model= train_model(X=x, Y=y)
with open(trained_model_file, 'wb+') as of:
    pickle.dump(trained_model, of)
