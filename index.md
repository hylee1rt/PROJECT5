For this project, we will be exploring dimensionality reduction methods and clustering techniques on high dimensional microbiome datasets. The first dataset is an adaptation of the data from the 2010 Fierer et al. paper, which can be found at (https://www.pnas.org/content/107/14/6477).

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


