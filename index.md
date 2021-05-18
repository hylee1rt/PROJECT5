# Phylum-level Clustering Analysis on Skin Microbial Communities

For this project, we will be exploring dimensionality reduction methods and clustering techniques on high dimensional microbiome communities. The dataset we will work with is an adaptation of the data from the 2010 Fierer et al. paper, which can be found at (https://www.pnas.org/content/107/14/6477). 

First, we will import and preprocess our data. 

```python
otus_data = pd.read_excel('otus.xlsx')
otus_data = otus_data.drop(columns=['Taxonomy','Genus','#OTU ID'])
otus_data = otus_data.groupby(["Phylum"]).sum()
otus_data = otus_data.T
otus_data = otus_data.reset_index()
otus_data[['Null','Label']] = otus_data['index'].str.split(".",expand=True)
otus_data = otus_data.drop(columns=['index','Null'])
otus_data[['Indiv']] = otus_data['Label'].str[:2]
otus_data[['Key']] = otus_data['Label'].str[2:-3:]
cols = [' Acidobacteria',' Actinobacteria',' Bacteroidetes',' Cyanobacteria',' Firmicutes',' Fusobacteria',' Gemmatimonadetes',' Proteobacteria','Indiv','Key']
otus_data = otus_data[cols]
```

This will finally bring us to a table that looks like this:
![Untitled](https://user-images.githubusercontent.com/66886936/116767579-0e80aa80-a9ff-11eb-96f9-503d9677fd2b.jpg)


The heatmaps show that there correlations between some of our variables. 
![download (1)](https://user-images.githubusercontent.com/66886936/116767645-733c0500-a9ff-11eb-8cd6-838452edede8.png)
![download (3)](https://user-images.githubusercontent.com/66886936/116769500-acc63d80-aa0a-11eb-88c5-b3d50fb14639.png)


From the parallel coordinates plane, we can see that most of our data belongs to M2, M3, or M9, and there is a working pattern that we might be able to explore. 

For further analyses, we will filter out other individuals and keep M2, M3 and M9. We will also standardize the data and work with a table of relative frequencies. 
![download (2)](https://user-images.githubusercontent.com/66886936/116769290-f4e46080-aa08-11eb-9834-7f16e3316259.png)



### PCA

Next we will perform Principle Component Analysis (PCA) on the dataset to see if we could discover anything interesting. PCA is a dimensionality-reduction method that is often used to reduce the dimensionality of large data sets, by transforming a large set of variables into a smaller one that still contains most of the information in the large set. The goal of PCA is basically to find an angle where you can project those points down onto a flat piece of paper, in such a way that the variance (or spread) of the data is maximized.

The intuition behind this is that, if we think we might have "clusters" or "groups" in a high-dimensional dataset, we might be able to find a particular angle of view in the high dimensional space, ~smoosh~ the data down onto a 2-d surface and be able to see the clustering or grouping in the data. This can help us decide whether a particular type of model might be useful.


```python
plt.figure(figsize = [8,8])

colors = ['magenta','cyan','red']
options = ['M2', 'M3', 'M9']

for x in range(len(options)):
    idx = (otus_Y['Indiv']==options[x])
    plt.scatter(pca_otus[idx,0],pca_otus[idx,1],
                edgecolor='k', color = colors[int(x)],alpha=0.5)

m_patch = mpatches.Patch(color='magenta', alpha = 0.6, label='M2')
c_patch = mpatches.Patch(color='cyan', alpha = 0.6, label='M3')
r_patch = mpatches.Patch(color='red', alpha = 0.6, label='M9')

plt.xlabel('PC1: '  + str(int(round(pca.explained_variance_ratio_[0]*100))) + '%',fontsize=14)
plt.ylabel('PC2: ' + str(int(round(pca.explained_variance_ratio_[1]*100))) + '%',fontsize=14)
plt.legend(handles=[m_patch, c_patch, r_patch], labels = options)
plt.title("PCA on Bacterial Phyla Abundances by Individual", fontsize=12)
plt.savefig('otus_pca.png', bbox_inches = 'tight')
plt.show()
```

![download (5)](https://user-images.githubusercontent.com/66886936/116770270-5d830b80-aa10-11eb-9e5a-a2f14d713b90.png)

Another thing to look at is how much variance is explained by each Principal Component. The first two principle components seem to explain over 98% of our data.  

```python
pca.explained_variance_ratio_

array([9.59566288e-01, 2.93070077e-02, 5.99185498e-03, 4.99038685e-03,
       1.36521106e-04, 7.14345547e-06, 7.97462217e-07, 7.90823233e-33])
```
![download (6)](https://user-images.githubusercontent.com/66886936/116770424-e189c300-aa11-11eb-83c9-c54c2cab407a.png)

We will take a look at a couple machine learning techniques that can be used to test how well those methods detect and identify the significantly different groups in the dataset as well as predict which class an instance of the data belongs to. In our case, we will see how well the models identify the three individuals using the relative frequency of each of the microbiome communities collected. We will also explore if and how accurately a model can classify a data point to the right individual.

## Clustering 

### Agglomerative Heirarchical Clustering

AHC works by starting off with each data point in its own cluster. Beginning with a distance of zero, we gradually increase the distance, and when two clusters (which may only contain one point initially) fall within that distance of each other, they are merged into a single cluster. The "distance between clusters" depends on the type of linkage being used. This process continues until all the data has been merged into a single cluster. The result of this type of clustering can be summarized using dendrograms, which can then be inspected to determine the optimal number of clusters. 

We can try different types of linkages to compare, and it looks like our data might have two or four clusters from our dendograms. We can compare the dendograms and the accuracies of our four different linkage types. "Complete" and "ward" linkage types had the highest accuracy, and if we take a look at the confusion matrices, we can see that the "single" and "average" linkages were only able to detect two clusters. They grouped the points with the label "2" with all the "1"s. In other words, they were not able to detect the third individual labeled "2". 

|Linkage |     Dendogram        |Confusion Matrix| Accuracy|
|:------:|:--------------------:|:-----:|:-----:|
|Single  |![download (6)](https://user-images.githubusercontent.com/66886936/118550332-8dc5ec00-b72a-11eb-91c0-a4d3ae49c0aa.png)|<img width="85" alt="single" src="https://user-images.githubusercontent.com/66886936/118563801-3c732800-b73d-11eb-8ba0-561e21e1b4f0.png">|0.68|
|Average |![download (7)](https://user-images.githubusercontent.com/66886936/118550377-99191780-b72a-11eb-8d2d-ce476c6036b5.png)|<img width="85" alt="average" src="https://user-images.githubusercontent.com/66886936/118564059-a25faf80-b73d-11eb-9556-7f723cce8e9c.png">|0.67|
|Complete |![download (8)](https://user-images.githubusercontent.com/66886936/118550420-a7673380-b72a-11eb-98a6-8dc58734a6c2.png)|<img width="93" alt="complete" src="https://user-images.githubusercontent.com/66886936/118564737-d9829080-b73e-11eb-9964-90f750f3b7b6.png">|0.75|
|Ward     |![download (9)](https://user-images.githubusercontent.com/66886936/118550460-b4842280-b72a-11eb-9ea8-df1c5f1ba0d1.png)|<img width="91" alt="ward" src="https://user-images.githubusercontent.com/66886936/118564831-059e1180-b73f-11eb-9285-b1f03877f4d7.png">|0.75|


However, since we already know how many clusters there are (there are 3 individuals in our data), we will specify the number of clusters. Using 3 as the number of clusters and the linkage type "complete," we can achieve a 75% accuracy with this clustering method. This means that our model was able to correctly group the data into three clusters that belong to three different individuals 75% of the time. We can also see the incorrectly clustered data points from the confusion matrices. 

```python
from sklearn.cluster import AgglomerativeClustering as AC
from sklearn.metrics import accuracy_score

n = 3
ac = AC(n_clusters=n,linkage='complete')
ac_cl = ac.fit_predict(X)
accuracy = accuracy_score(y,ac_cl)
print(accuracy)
0.75
```
The TSNE plot below visualizes the clusters that our model (n_clusters = 3 and linkage = "complete") predicted, and we can see that there are a couple of green points where all purple points are, as we can also see from the confusion matrix that the model labeled one two of the 0's as 1 and 2. 

![download](https://user-images.githubusercontent.com/66886936/118583999-fb423e80-b763-11eb-93d5-a76cf5e99582.png)


## Classification

### k-nearest neighbors (KNN)

Classification refers to types of supervised learning. With supervised learning, the model is being provided with information about what the correct answers are during the model training process. Since our data has the labels for each individual (we know who each record of microbiome communities belongs to), we can explore how classification methods perform on our dataset. KNN classification works by runs a mathematical formula to compute the distance between each data point and the test data, selecting the specified number points (K) closest to it, then finds the probability of these points being similar to the test data and classifies it based on which points share the highest probabilities.

For this method, we will split our data into training and test sets to compare the predictions to the actual classes. 

```python
from sklearn.neighbors import KNeighborsClassifier as KNN
from sklearn.model_selection import train_test_split as tts

def CompareClasses(actual, predicted, names=None):
    '''Function returns a confusion matrix, and overall accuracy given:
            Input:  actual - a list of actual classifications
                    predicted - a list of predicted classifications
                    names (optional) - a list of class names
    '''
    
    import pandas as pd
    accuracy = sum(actual==predicted)/actual.shape[0]
    classes = pd.DataFrame(columns=['Actual','Predicted'])
    classes['Actual'] = actual
    classes['Predicted'] = predicted
    conf_mat = pd.crosstab(classes['Predicted'],classes['Actual'])
    # Relabel the rows/columns if names was provided
    if type(names) != type(None):
        conf_mat.index=y_names
        conf_mat.index.name='Predicted'
        conf_mat.columns=y_names
        conf_mat.columns.name = 'Actual'
    print('Accuracy = ' + format(accuracy, '.2f'))
    return conf_mat, accuracy
    
    
Xtrain,Xtest,ytrain,ytest = tts(X, y, test_size=0.25, random_state=2021)

knn = KNN(n_neighbors=12)
knn.fit(Xtrain,ytrain)
knn_cl = knn.predict(Xtest)
CompareClasses(ytest,knn_cl)

Accuracy = 0.88
(Actual      0  1  2
 Predicted          
 0          13  0  0
 1           0  4  1
 2           0  2  6, 0.8846153846153846)
```

Here, we have an accuracy of 88%. The resulting confusion matrix shows that some 1's were classified as 2, and some 2's were classified as 1, but the model was able to classify all 0's correctly. The pca plot below of our test data shows that group 1 and 2 are closer together than group 0, which can be explained by how our models had trouble distinguishing labels 1 and 2. 


![download (1)](https://user-images.githubusercontent.com/66886936/118585288-920ffa80-b766-11eb-828b-b3837e92bd86.png)







