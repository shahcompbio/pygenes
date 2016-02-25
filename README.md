# PyGenes

PyGenes is an in memory searchable gene database for python, implemented in c++ for speed.  The current implementation is primarily aimed at providing an interface for ensembl gtf files.

## Usage

Create a pygenes database and load from an ensembl gtf file:

```
gene_models = pygenes.GeneModels()
gene_models.load_ensembl_gtf('tests/Homo_sapiens.NCBI36.54.test.gtf')
``` 

The database can be serialized to a binary format:

```
gene_models.save_binary('gene_models.binary')

gene_models2 = pygenes.GeneModels()
gene_models2.load_binary('gene_models.binary')
```

Get a gene by gene id:

```
gene = gene_models.get_gene('ENSG00000101596')
```

Find either overlapping or contained genes within a region:

```
gene_models.find_overlapping_genes('18', 2700000, 2800000)
gene_models.find_contained_genes('18', 2700000, 2800000)
```

Calculate the location in the gene as one of 'utr3p', 'coding', 'intron', 'utr5p', 'upstream', 'downstream', 'utr':

```
gene_models.calculate_gene_location('ENSG00000101596', 2792681)
```

Calculate the location in the genome of a position in a transcript:

```
gene_models.calculate_genomic_position('ENST00000320876', 461)
```

Calculate the regions in the genome of a region in a transcript:

```
gene_models.calculate_genomic_regions('ENST00000320876', 461, 796)
```

