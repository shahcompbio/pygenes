# distutils: language = c++
# cython: language_level=3, c_string_type=unicode, c_string_encoding=utf8


from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "genemodels.cpp":
    cdef cppclass CRegion:
        int start
        int end
    cdef cppclass CGene:
        string id
        string name
        string source
        string chromosome
        string strand
        int start
        int end
    cdef cppclass CGeneModels:
        void CGeneModels() except +
        void LoadEnsemblGTF(string gtf_filename) except +
        CGene GetGene(string geneID) except +
        string GetTranscriptGene(string transcriptID) except +
        void FindOverlappingGenes(string chromosome, int start, int end, vector[string] genes) except +
        void FindContainedGenes(string chromosome, int start, int end, vector[string] genes) except +
        void FindNearestGenes(string chromosome, int position, vector[string] genes) except +
        string CalculateGeneLocation(string geneID, int position) except +
        int CalculateGenomicPosition(string transcriptID, int position) except +
        void CalculateGenomicRegions(string transcriptID, int start, int end, vector[CRegion]& regions) except +


cdef extern from "IntervalTree/IntervalTree.h":
    cdef cppclass CInterval[T]:
        int start
        int stop
        T value
    cdef cppclass CIntervalTree[T]:
        void CIntervalTree(vector[CInterval[T]] intervals) except +
        void FindOverlapping(int start, int stop, vector[T] overlapping) except +
        void FindContained(int start, int stop, vector[T] contained) except +
        void FindNearest(int position, vector[T] nearest) except +


class Region:
    def __init__(self, start, end):
        self.start = start
        self.end = end


class Gene:
    def __init__(self, id, name, source, chromosome, strand, start, end):
        self.id = id
        self.name = name
        self.source = source
        self.chromosome = chromosome
        self.strand = strand
        self.start = start
        self.end = end


cdef class IntervalTree:
    cdef CIntervalTree[int] *c_interval_tree

    def __cinit__(self, intervals):
        cdef CInterval[int] c_interval
        cdef vector[CInterval[int]] c_intervals = vector[CInterval[int]]()

        for interval in intervals:
            if len(interval) != 3:
                raise ValueError('excpected tuple of size 3')

            c_interval.value = interval[0]
            c_interval.start = interval[1]
            c_interval.stop = interval[2]

            c_intervals.push_back(c_interval)

        self.c_interval_tree = new CIntervalTree[int](c_intervals)

    def __dealloc__(self):
        del self.c_interval_tree

    def find_overlapping(self, start, stop):
        cdef vector[int] overlapping
        self.c_interval_tree.FindOverlapping(start, stop, overlapping)
        return overlapping

    def find_contained(self, start, stop):
        cdef vector[int] contained
        self.c_interval_tree.FindContained(start, stop, contained)
        return contained

    def find_nearest(self, position):
        cdef vector[int] nearest
        self.c_interval_tree.FindNearest(position, nearest)
        return nearest


cdef class GeneModels:
    cdef CGeneModels *c_gene_models

    def __cinit__(self):
        self.c_gene_models = new CGeneModels()

    def __dealloc__(self):
        del self.c_gene_models

    def load_ensembl_gtf(self, str gtf_filename):
        self.c_gene_models.LoadEnsemblGTF(gtf_filename)

    def get_gene(self, gene_id):
        cdef CGene gene = self.c_gene_models.GetGene(gene_id)
        return Gene(
            gene.id, gene.name, gene.source, gene.chromosome,
            gene.strand, gene.start, gene.end)

    def get_transcript_gene(self, transcript_id):
        return self.c_gene_models.GetTranscriptGene(transcript_id)

    def find_overlapping_genes(self, chromosome, start, end):
        cdef vector[string] genes
        self.c_gene_models.FindOverlappingGenes(chromosome, start, end, genes)
        return genes

    def find_contained_genes(self, chromosome, start, end):
        cdef vector[string] genes
        self.c_gene_models.FindContainedGenes(chromosome, start, end, genes)
        return genes

    def find_nearest_genes(self, chromosome, position):
        cdef vector[string] genes
        self.c_gene_models.FindNearestGenes(chromosome, position, genes)
        return genes

    def calculate_gene_location(self, gene_id, position):
        return self.c_gene_models.CalculateGeneLocation(gene_id, position)

    def calculate_genomic_position(self, transcript_id, position):
        return self.c_gene_models.CalculateGenomicPosition(transcript_id, position)

    def calculate_genomic_regions(self, transcript_id, start, end):
        cdef vector[CRegion] regions
        self.c_gene_models.CalculateGenomicRegions(transcript_id, start, end, regions)
        return [Region(r.start, r.end) for r in regions]
