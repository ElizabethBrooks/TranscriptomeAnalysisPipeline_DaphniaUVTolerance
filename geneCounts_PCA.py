import matplotlib as mpl
## agg backend is used to create plot as a .png file
mpl.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
#Retrieve track features and popularities and store in data frame
cols =  ['TrackID', 'Popularity', 'FeatureID', 'Danceability', 'Energy', 'Key', 'Loudness', 'Mode', 
         'Speechiness', 'Acousticness', 'Instrumentalness', 'Livenss', 
         'Valence', 'Tempo', 'TimeSignature', 'Duration', 'Bin']
data = pd.read_csv('trackFeatures_totalTrackFeaturesMapped_shuffled_binned.csv', names=cols)
y = data['Bin']          # Split off classifications
X = data.ix[:, 'Danceability':'Duration'] # Split off features
####
#Normalize data
X_norm = (X - X.min())/(X.max() - X.min())
####
#Principle component analysis
pca = sklearnPCA(n_components=2) #2-dimensional PCA
transformed = pd.DataFrame(pca.fit_transform(X_norm))
#Plot normalized data
plt.scatter(transformed[y==0][0], transformed[y==0][1], label='Popularity = 0', c='red')
plt.scatter(transformed[y==1][0], transformed[y==1][1], label='0 < Popularity < 25', c='green')
plt.scatter(transformed[y==2][0], transformed[y==2][1], label='24 < Popularity < 50', c='orange')
plt.scatter(transformed[y==3][0], transformed[y==3][1], label='49 < Popularity < 75', c='purple')
plt.scatter(transformed[y==4][0], transformed[y==4][1], label='74 < Popularity < 100', c='pink')
plt.scatter(transformed[y==5][0], transformed[y==5][1], label='Popularity = 100', c='blue')
#display plot
#plt.show()
#Save plot
plt.savefig('scatterPlot_PCA.png', bbox_inches='tight')