import os
import unittest
import pygenes


class pygenes_test(unittest.TestCase):
    
    def setUp(self):
        
        self.gtf_filename = os.path.join(os.path.dirname(__file__), 'Homo_sapiens.NCBI36.54.test.gtf')
    
    def test_simple(self):
        
        gene_models = pygenes.GeneModels()
        gene_models.load_ensembl_gtf(self.gtf_filename)
        gene_models.save_binary('gene_models.binary')
        
        gene_models = pygenes.GeneModels()
        gene_models.load_binary('gene_models.binary')
        
        gene = gene_models.get_gene('ENSG00000101596')
        self.assertEqual(gene.id, 'ENSG00000101596')
        self.assertEqual(gene.name, 'SMCHD1')
        self.assertEqual(gene.source, 'protein_coding')
        self.assertEqual(gene.chromosome, '18')
        self.assertEqual(gene.strand, '+')
        self.assertEqual(gene.start, 2690857)
        self.assertEqual(gene.end, 2792925)
        
        self.assertEqual(gene_models.get_transcript_gene('ENST00000379020'), 'ENSG00000065665')
        
        self.assertEqual(set([a for a in gene_models.find_overlapping_genes('18', 2700000, 2800000)]),
                         set(['ENSG00000101596', 'ENSG00000209536', 'ENSG00000207034']))
        
        self.assertEqual(set([a for a in gene_models.find_contained_genes('18', 2700000, 2800000)]),
                         set(['ENSG00000209536', 'ENSG00000207034']))
        
        self.assertEqual(set([a for a in gene_models.find_nearest_genes('18', 2810000)]),
                         set(['ENSG00000101596']))

        self.assertEqual(gene_models.calculate_gene_location('ENSG00000101596', 2792681), 'utr3p')
        self.assertEqual(gene_models.calculate_gene_location('ENSG00000101596', 2767810), 'coding')
        self.assertEqual(gene_models.calculate_gene_location('ENSG00000101596', 2767355), 'intron')
        self.assertEqual(gene_models.calculate_gene_location('ENSG00000101596', 2693692), 'utr5p')
        self.assertEqual(gene_models.calculate_gene_location('ENSG00000101596', 2689083), 'upstream')
        self.assertEqual(gene_models.calculate_gene_location('ENSG00000101596', 2805042), 'downstream')
        self.assertEqual(gene_models.calculate_gene_location('ENSG00000209536', 2764898), 'utr')
        
        self.assertEqual(gene_models.calculate_genomic_position('ENST00000320876', 461), 2656878)
        
        self.assertEqual([(a.start, a.end) for a in gene_models.calculate_genomic_regions('ENST00000320876', 461, 796)],
                         [(2656878, 2657030), (2663280, 2663362), (2664014, 2664113)])
        
        
if __name__ == '__main__':
    unittest.main()
    

