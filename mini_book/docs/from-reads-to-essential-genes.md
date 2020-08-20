---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.8'
    jupytext_version: 1.4.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# From reads and insertions to essentiality prediction using a Regression Model

```{code-cell} ipython3
import pandas as pd
import pylab as pl
import numpy as np
import seaborn as sns
import scipy.optimize as opt
from sklearn import preprocessing
%matplotlib inline 
import matplotlib.pyplot as plt
```

```{code-cell} ipython3
import os
script_dir = os.path.dirname('__file__') #<-- absolute dir the script is in
rel_path_data_insertions="Python_scripts/Data_Files/wt-truncated-insertions-per-gene.xlsx"
rel_path_data_reads="Python_scripts/Data_Files/WT-benoit-truncated-reads-per-gene.xlsx"


abs_path_data_reads = os.path.join(script_dir, rel_path_data_reads)
abs_path_data_insertions = os.path.join(script_dir, rel_path_data_insertions)

# os.chdir('../') #<-- for binder os.chdir('../')
```

```{code-cell} ipython3
data_reads=pd.read_excel(abs_path_data_reads)
data_insertions=pd.read_excel(abs_path_data_insertions)
```

```{code-cell} ipython3
data_reads_pd=data_reads.iloc[:,0:3]
data_reads_pd.columns=['Gene_name','Essentiality','reads-truncated']
data_reads_pd.head()
```

```{code-cell} ipython3
data_insertions_pd=data_insertions.iloc[:,0:3]
data_insertions_pd.columns=['Gene_name','Essentiality','insertions-truncated']
data_insertions_pd.head()
```

```{code-cell} ipython3
data=data_reads_pd.copy()
data['insertions']=data_insertions_pd['insertions-truncated']
```

```{code-cell} ipython3
data
```

## Visualising how the input data is distributed regarding the type of essentiality of the genes

```{code-cell} ipython3
fig, axes=plt.subplots(1,2)
plt.subplots_adjust(wspace=0.5)
sns.boxplot(x='Essentiality',y='insertions',data=data,ax=axes[0])
# sns.boxplot(x='Essentiality',y='insertions-non-truncated',data=data,ax=axes[0])
axes[0].set_ylim([0,100])
sns.boxplot(x='Essentiality',y='reads-truncated',data=data,ax=axes[1])
axes[1].set_ylim([0,300])
```

## Remove outliers 

```{code-cell} ipython3
sns.scatterplot(y='insertions',x=np.arange(0,len(data)),data=data,hue='Essentiality')
# plt.ylim([0,20])
plt.figure()
sns.scatterplot(y='reads-truncated',x=np.arange(0,len(data)),data=data,hue='Essentiality')
# plt.ylim([0,260])
```

### Make use of the IQR to remove the outliers

- IQR is part of Descriptive statistics and also called as midspead , middle 50%

- IQR is first Quartile minus the Third Quartile (Q3-Q1)

```{code-cell} ipython3
y_insert=data['insertions']
removed_outliers_insertions = y_insert.between(y_insert.quantile(.05), y_insert.quantile(.95)) # any value bellow  Q1-0.05*IQR or above Q3+0.95*IQR is an outlier

y_reads=data['reads-truncated']
removed_outliers_reads = y_reads.between(y_reads.quantile(.05), y_reads.quantile(.95))
```

```{code-cell} ipython3
:tags: []

fig, axes=plt.subplots(1,2)
plt.subplots_adjust(wspace=0.5,right=1.2)

axes[0].plot(y_insert,alpha=0.3)
axes[0].plot(y_insert[removed_outliers_insertions],label='removed-outliers')
axes[0].legend()
axes[0].set_ylabel('Insertions after truncation')
axes[1].plot(y_reads,alpha=0.3)
axes[1].plot(y_reads[removed_outliers_reads],label='removed-outliers')
axes[1].legend()
axes[1].set_ylabel('Reads after truncation')
```

```{code-cell} ipython3
data['insertions-non-outliers']=y_insert[removed_outliers_insertions]
data['reads-non-outliers']=y_reads[removed_outliers_reads]
```

```{code-cell} ipython3
fig, axes=plt.subplots(1,4)
plt.subplots_adjust(wspace=2.6,right=1.5)

sns.boxplot(x='Essentiality',y='insertions-non-outliers',data=data,ax=axes[0])
axes[0].set_title('With the outliers removed')
sns.boxplot(x='Essentiality',y='insertions',data=data,ax=axes[1])
axes[1].set_title('Original')
sns.boxplot(x='Essentiality',y='reads-non-outliers',data=data,ax=axes[2])
axes[2].set_title('With the outliers removed')
sns.boxplot(x='Essentiality',y='reads-truncated',data=data,ax=axes[3])
axes[3].set_title('original')
```

```{code-cell} ipython3
data_remove_outliers=data.copy()
# data_remove_outliers=data_remove_outliers.drop(columns=['reads-non-truncated','insertions-non-truncated'])
data_remove_outliers.fillna(0,inplace=True)
```

```{code-cell} ipython3
data_remove_outliers.head()
```

```{code-cell} ipython3
data_remove_outliers.groupby('Essentiality').describe()
```

```{code-cell} ipython3
data_remove_outliers[data_remove_outliers['Gene_name']=='CDC24']
```

```{code-cell} ipython3
sns.pairplot(data=data_remove_outliers,hue='Essentiality',vars=['reads-non-outliers','insertions-non-outliers'])
```

```{code-cell} ipython3
X = np.asarray(data_remove_outliers[['insertions-non-outliers', 'reads-non-outliers']])
X[0:5]
```

```{code-cell} ipython3
y = np.asarray(data['Essentiality'])
y [0:30]
```

```{code-cell} ipython3
from sklearn import preprocessing
X = preprocessing.StandardScaler().fit(X).transform(X)
X[0:5]
```

```{code-cell} ipython3
:tags: []

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split( X, y, test_size=0.1, random_state=4)
print ('Train set:', X_train.shape,  y_train.shape)
print ('Test set:', X_test.shape,  y_test.shape)
```

```{code-cell} ipython3
X_train[np.isnan(X_train)]
```

```{code-cell} ipython3
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
LR = LogisticRegression(C=0.01, solver='liblinear').fit(X_train,y_train)
LR
```

```{code-cell} ipython3
yhat = LR.predict(X_test)
yhat[0:4]
```

```{code-cell} ipython3
:tags: []

from sklearn.metrics import classification_report, confusion_matrix
import itertools
def plot_confusion_matrix(cm, classes,
                          normalize=False,
                          title='Confusion matrix',
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
print(confusion_matrix(y_test, yhat, labels=[1,0]))
```

```{code-cell} ipython3
:tags: []

cnf_matrix = confusion_matrix(y_test, yhat, labels=[1,0])
np.set_printoptions(precision=2)


# Plot non-normalized confusion matrix
plt.figure()
plot_confusion_matrix(cnf_matrix, classes=['Essential','Non essential'],normalize= False,  title='Confusion matrix')
```

```{code-cell} ipython3
:tags: []

print (classification_report(y_test, yhat))
```

### Predicting the probabilities for each gene to be essential or not

```{code-cell} ipython3
yprob=LR.predict_proba(X)# probability of cdc24 to be essential according this model
yprob[35]
```

```{code-cell} ipython3
data_remove_outliers['probability of being essential']=yprob[:,1]
```

```{code-cell} ipython3
sns.distplot(data_remove_outliers[data_remove_outliers['Essentiality']==0]['probability of being essential'],label='non-essentials',norm_hist=True)
sns.distplot(data_remove_outliers[data_remove_outliers['Essentiality']==1]['probability of being essential'],label='true-essentials',norm_hist=True)
plt.vlines(x=0.3,ymin=0,ymax=60,linewidth=2,alpha=0.5)
plt.vlines(x=0.4,ymin=0,ymax=60,linewidth=2,alpha=0.5)
plt.hlines(y=60,xmin=0.3,xmax=0.4,linewidth=2,alpha=0.5)
plt.fill_between(x=[0.3, 0.4],y1=60,color='black',alpha=0.2,label='76% of changes to be essential')
plt.legend()
```

```{code-cell} ipython3
ratio_essentials=[]
ratio_nonessentials=[]
prob=np.arange(0,1,0.01)
for i in prob:
    ratio_essentials.append(np.sum(y[np.where(yprob[:,1]>i)[0]])/len(yprob[:,1]>i)) # the sum symbolizes the number of 1's which are the true essential proteins 
    ratio_nonessentials.append(1-np.sum(y[np.where(yprob[:,0]>i)[0]])/len(y[np.where(yprob[:,0]>i)[0]]))

plt.scatter(prob,ratio_nonessentials,alpha=0.3,label= 'Class 0 Probab')
plt.scatter(prob,ratio_essentials,alpha=0.3,label='Class I Probab')
plt.xlabel('Probabilities')
plt.ylabel('ratio of true classification')
plt.legend()
```

```{code-cell} ipython3
index_prone2essential=[]
bound=0.3
for i in np.arange(0,len(yprob)):
    if yprob[i,1] > bound : 
        index_prone2essential.append(i)
ratio_true_essentials=np.sum(data_remove_outliers.iloc[index_prone2essential,:]['Essentiality'])/data_remove_outliers['Essentiality'].value_counts()[1]
```

### Summary from the model

```{code-cell} ipython3
:tags: []

print('- The essential genes represents a ',data_remove_outliers['Essentiality'].value_counts()[1]/len(data_remove_outliers) * 100,'%', 'of the population of genes.')

print('- The ratio of true essentials contained in the genes that have more than',bound,'probability of being essential is =', ratio_true_essentials *100,'%.','This represents',(len(index_prone2essential))/len(data_remove_outliers) * 100,'%', 'of all genes. In other words, with this regression model we can trust that a gene has', ratio_true_essentials *100,'%' ,'changes of being essential if the probability given by the model is higher than',bound)


print('- For the genes that has less than ',bound,'Pr to be essential , then I can assure they are non essential, which are ', len(data_remove_outliers)-len(index_prone2essential),'genes. That represents the',(len(data_remove_outliers)-len(index_prone2essential))/data_remove_outliers['Essentiality'].value_counts()[0] * 100,'%','of the total number of non essential genes.')
```

```{code-cell} ipython3
plt.hist(y[np.where(yprob[:,1]<0.3)]);
```

```{code-cell} ipython3
sns.distplot(yprob[:,0])
sns.distplot(yprob[:,1])
```

```{code-cell} ipython3

```
