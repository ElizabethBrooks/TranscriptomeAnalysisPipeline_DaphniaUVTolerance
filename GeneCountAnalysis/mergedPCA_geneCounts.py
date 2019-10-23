import matplotlib as mpl
## agg backend is used to create plot as a .png file
mpl.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA as sklearnPCA
#Specify header tags for csv file of gene counts
cols =  ['gene', 'Y05_VIS_Pool1_T', 'Y05_VIS_Pool2_T', 'Y05_VIS_Pool3_T', 'Y05_UV_Pool1_T', 
'Y05_UV_Pool2_T', 'Y05_UV_Pool3_T', 'Y023_5_VIS_Pool1_T', 'Y023_5_VIS_Pool2_T', 'Y023_5_VIS_Pool3_T', 
'Y023_5_UV_Pool1_T', 'Y023_5_UV_Pool2_T', 'Y023_5_UV_Pool3_T', 'E05_VIS_Pool1_T', 'E05_VIS_Pool2_T', 
'E05_VIS_Pool3_T', 'E05_UV_Pool1_T', 'E05_UV_Pool2_T', 'E05_UV_Pool3_T', 'R2_VIS_Pool1_T', 
'R2_VIS_Pool2_T', 'R2_VIS_Pool3_T', 'R2_UV_Pool1_T', 'R2_UV_Pool2_T', 'R2_UV_Pool3_T', 
'Y05_VIS_Pool1_H', 'Y05_VIS_Pool2_H', 'Y05_VIS_Pool3_H', 'Y05_UV_Pool1_H', 'Y05_UV_Pool2_H', 
'Y05_UV_Pool3_H', 'Y023_5_VIS_Pool1_H', 'Y023_5_VIS_Pool2_H', 'Y023_5_VIS_Pool3_H', 
'Y023_5_UV_Pool1_H', 'Y023_5_UV_Pool2_H', 'Y023_5_UV_Pool3_H', 'E05_VIS_Pool1_H', 'E05_VIS_Pool2_H', 
'E05_VIS_Pool3_H', 'E05_UV_Pool1_H', 'E05_UV_Pool2_H', 'E05_UV_Pool3_H', 'R2_VIS_Pool1_H', 
'R2_VIS_Pool2_H', 'R2_VIS_Pool3_H', 'R2_UV_Pool1_H', 'R2_UV_Pool2_H', 'R2_UV_Pool3_H']
#Retrieve gene counts and store in data frame
data = pd.read_csv('../../final_merged_counts.csv', names=cols)
#Specify target sample tags
samples =  ['Y05_VIS_Pool1_T', 'Y05_VIS_Pool2_T', 'Y05_VIS_Pool3_T', 'Y05_UV_Pool1_T', 
'Y05_UV_Pool2_T', 'Y05_UV_Pool3_T', 'Y023_5_VIS_Pool1_T', 'Y023_5_VIS_Pool2_T', 'Y023_5_VIS_Pool3_T', 
'Y023_5_UV_Pool1_T', 'Y023_5_UV_Pool2_T', 'Y023_5_UV_Pool3_T', 'E05_VIS_Pool1_T', 'E05_VIS_Pool2_T', 
'E05_VIS_Pool3_T', 'E05_UV_Pool1_T', 'E05_UV_Pool2_T', 'E05_UV_Pool3_T', 'R2_VIS_Pool1_T', 
'R2_VIS_Pool2_T', 'R2_VIS_Pool3_T', 'R2_UV_Pool1_T', 'R2_UV_Pool2_T', 'R2_UV_Pool3_T', 
'Y05_VIS_Pool1_H', 'Y05_VIS_Pool2_H', 'Y05_VIS_Pool3_H', 'Y05_UV_Pool1_H', 'Y05_UV_Pool2_H', 
'Y05_UV_Pool3_H', 'Y023_5_VIS_Pool1_H', 'Y023_5_VIS_Pool2_H', 'Y023_5_VIS_Pool3_H', 
'Y023_5_UV_Pool1_H', 'Y023_5_UV_Pool2_H', 'Y023_5_UV_Pool3_H', 'E05_VIS_Pool1_H', 'E05_VIS_Pool2_H', 
'E05_VIS_Pool3_H', 'E05_UV_Pool1_H', 'E05_UV_Pool2_H', 'E05_UV_Pool3_H', 'R2_VIS_Pool1_H', 
'R2_VIS_Pool2_H', 'R2_VIS_Pool3_H', 'R2_UV_Pool1_H', 'R2_UV_Pool2_H', 'R2_UV_Pool3_H']
#Split off targets
y = data.loc[;, samples]
#Split off features
X = data.loc[;, ['gene']]
####
#Standardizing the features
X_stand = StandardScaler().fit_transform(X)
#Normalize the features
#X_norm = (X - X.min())/(X.max() - X.min())
####
#Principle component analysis
pca = sklearnPCA(n_components=2) #2-dimensional PCA
#transformed = pd.DataFrame(pca.fit_transform(X_norm))
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data=principalComponents, columns=['PC 1', 'PC 2'])
finalDf = pd.concat([principalDf, df[samples]], axis=1)
#Plot data
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('PC 1', fontsize=15)
ax.set_ylabel('PC 2', fontsize=15)
ax.set_title('Merged Gene Count PCA', fontsize=20)
targets = ['_T', '_H']
colors = ['r', 'b']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['target'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'PC 1']
               , finalDf.loc[indicesToKeep, 'PC 2']
               , c=color
               , s=50)
ax.legend(targets)
ax.grid()
#display plot
#fig.show()
#Save plot
fig.savefig('geneCounts_merged_PCA.png', bbox_inches='tight')