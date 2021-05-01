For this project, we will be exploring dimensionality reduction methods and clustering techniques on high dimensional microbiome communities and tumor biopsy datasets. The first dataset we will work with is an adaptation of the data from the 2010 Fierer et al. paper, which can be found at (https://www.pnas.org/content/107/14/6477).

## Phylum-level Clustering Analysis on Skin Microbial Communities

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


From the parallel coordinates plane, we can see that most of our data belongs to M2, M3, or M9, and there is a working pattern that we might be able to explore. We will filter out other individuals and keep M2, M3 and M9. We will also standardize the data and work with a table of relative frequencies. 
![download (2)](https://user-images.githubusercontent.com/66886936/116769290-f4e46080-aa08-11eb-9834-7f16e3316259.png)



### PCA

Next we will perform Principle Component Analysis (PCA) on the dataset to see if we could discover anything interesting. 

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

![download (4)](https://user-images.githubusercontent.com/66886936/116769723-49d5a600-aa0c-11eb-9734-84387557dab6.png)
