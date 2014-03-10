
#include "IntervalTree/IntervalTree.h"

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/pure_virtual.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/sum.hpp>

#include <fstream>
#include <vector>

#include "boost/serialization/unordered_map.hpp"
#include "boost/serialization/unordered_set.hpp"

using namespace std;
using namespace boost;


template<typename T>
void Print(const vector<T>& intervals)
{
	for (typename vector<T>::const_iterator intervalIter = intervals.begin(); intervalIter != intervals.end(); intervalIter++)
	{
		cout << *intervalIter << endl;
	}	
}


template <typename T>
inline python::list convert_vector(const vector<T>& v)
{
	python::list converted;
	for (typename vector<T>::const_iterator iter = v.begin(); iter != v.end(); iter++)
	{
		converted.append(*iter);
	}
	return converted;
}


struct Region
{
	Region() : start(-1), end(-1) {}
	Region(int start, int end) : start(start), end(end) {}
	
	int start;
	int end;
	
	int GetStart() const { return start; }
	int GetEnd() const { return end; }
	int GetLength() const { return end - start + 1; }
	
	bool operator==(const Region& other) const
	{
		return start == other.start && end == other.end;
	}
	
	bool operator<(const Region& other) const
	{
		return start < other.start;
	}
	
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & start;
		ar & end;
	}
};


struct GenomicRegion : public Region
{
	string chromosome;
	string strand;
	
	string GetChromosome() const { return chromosome; }
	string GetStrand() const { return strand; }
	
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & chromosome;
		ar & strand;
        ar & serialization::base_object<Region>(*this);
	}
};


struct Gene : public GenomicRegion
{
	string id;
	string name;
	string source;
	
	string GetID() const { return id; }
	string GetName() const { return name; }
	string GetSource() const { return source; }
	
	bool operator==(const Gene& other) const
	{
		return id == other.id;
	}
	
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & id;
		ar & name;
		ar & source;
        ar & serialization::base_object<GenomicRegion>(*this);
	}
};


class GeneModels
{
public:
	void LoadEnsemblGTF(const string& gtfFilename)
	{
		unordered_map<string, accumulators::accumulator_set<int, accumulators::stats<accumulators::tag::min, accumulators::tag::max > > > geneRegionAcc;
		
		// Open clusters file
		ifstream gtfFile(gtfFilename.c_str());
		if (!gtfFile)
		{
			stringstream errorStr;
			errorStr << "File " << gtfFilename << " not found";
			throw std::invalid_argument(errorStr.str());
		}
		
		// Parse file contents
		string line;
		int lineNumber = 0;
		while (getline(gtfFile, line))
		{
			lineNumber++;
			
			if (line.substr(0, 2) == "#!")
			{
				continue;
			}
			
			if (line.length() == 0)
			{
				stringstream errorStr;
				errorStr << "Empty gtf line " << lineNumber << " of " << gtfFilename;
				throw std::invalid_argument(errorStr.str());
			}
			
			vector<string> gtfFields;
			split(gtfFields, line, is_any_of("\t"));
			
			if (gtfFields.size() < 9)
			{
				stringstream errorStr;
				errorStr << "Error: Format error for gtf line " << lineNumber << " of " << gtfFilename;
				throw std::invalid_argument(errorStr.str());
			}
			
			const string& chromosome = gtfFields[0];
			const string& source = gtfFields[1];
			const string& featureType = gtfFields[2];
			int start = lexical_cast<int>(gtfFields[3]);
			int end = lexical_cast<int>(gtfFields[4]);
			const string& strand = gtfFields[6];
			
			vector<string> featureFields;
			split(featureFields, gtfFields[8], is_any_of(";"));
			
			string geneID;
			string transcriptID;
			string geneName;
			int exonNumber = -1;
			for (vector<string>::iterator featureIter = featureFields.begin(); featureIter != featureFields.end(); featureIter++)
			{
				trim(*featureIter);
				
				if (featureIter->empty())
				{
					continue;
				}
				
				int firstSpace = featureIter->find_first_of(' ');
				
				string key = featureIter->substr(0, firstSpace);
				string value = featureIter->substr(firstSpace + 2, featureIter->size() - firstSpace - 3);			
				
				if (key == "gene_id")
				{
					geneID = value;
				}
				else if (key == "transcript_id")
				{
					transcriptID = value;
				}
				else if (key == "gene_name")
				{
					geneName = value;
				}
				else if (key == "exon_number")
				{
					exonNumber = lexical_cast<int>(value);
				}
				
				assert(geneID.empty());
				assert(transcriptID.empty());
				assert(geneName.empty());
				assert(exonNumber != -1);
			}
			
			mGenes[geneID].id = geneID;
			mGenes[geneID].name = geneName;
			mGenes[geneID].source = source;
			mGenes[geneID].chromosome = chromosome;
			mGenes[geneID].strand = strand;
			geneRegionAcc[geneID](start);
			geneRegionAcc[geneID](end);
			
			mGeneTranscripts[geneID].insert(transcriptID);
			mTranscriptGene[transcriptID] = geneID;
			
			if (featureType == "exon")
			{
				mExons[transcriptID].push_back(Region(start, end));
			}
			else if (featureType == "CDS")
			{
				mCDSs[transcriptID].push_back(Region(start, end));
			}
			else if (featureType == "start_codon")
			{
				mStartCodon[transcriptID] = Region(start, end);
			}
			else if (featureType == "stop_codon")
			{
				mStopCodon[transcriptID] = Region(start, end);
			}
		}
		
		unordered_map<string,vector<Interval<string> > > geneIntervals;
		
		for (unordered_map<string,Gene>::const_iterator geneIter = mGenes.begin(); geneIter != mGenes.end(); geneIter++)
		{
			const string geneID = geneIter->first;
			
			int geneStart = accumulators::min(geneRegionAcc[geneID]);
			int geneEnd = accumulators::max(geneRegionAcc[geneID]);
			
			mGenes[geneID].start = geneStart;
			mGenes[geneID].end = geneEnd;
			
			geneIntervals[mGenes[geneID].chromosome].push_back(Interval<string>(geneStart, geneEnd, geneID));
		}
		
		mGeneIntervalTrees.clear();
		for (unordered_map<string,vector<Interval<string> > >::iterator chromosomeIter = geneIntervals.begin(); chromosomeIter != geneIntervals.end(); chromosomeIter++)
		{
			mGeneIntervalTrees[chromosomeIter->first] = IntervalTree<string>(chromosomeIter->second);
		}
		
		for (unordered_map<string,vector<Region> >::iterator exonsIter = mExons.begin(); exonsIter != mExons.end(); exonsIter++)
		{
			sort(exonsIter->second.begin(), exonsIter->second.end());
			
			int length = 0;
			for (vector<Region>::const_iterator exonIter = exonsIter->second.begin(); exonIter != exonsIter->second.end(); exonIter++)
			{
				length += exonIter->end - exonIter->start + 1;
			}
			mTranscriptLength[exonsIter->first] = length;
		}
		
		for (unordered_map<string,vector<Region> >::iterator cdssIter = mCDSs.begin(); cdssIter != mCDSs.end(); cdssIter++)
		{
			sort(cdssIter->second.begin(), cdssIter->second.end());
		}
	}
	
	void SaveBinary(const string& filename)
	{
		ofstream outputFile(filename.c_str());
		archive::binary_oarchive outputArchive(outputFile);
		outputArchive << *this;
	}
	
	void LoadBinary(const string& filename)
	{
		ifstream inputFile(filename.c_str());
		archive::binary_iarchive inputArchive(inputFile);
		inputArchive >> *this;
	}
	
	Gene GetGene(const string& geneID)
	{
		return mGenes[geneID];
	}
	
	string GetTranscriptGene(const string& transcriptID)
	{
		return mTranscriptGene[transcriptID];
	}
	
	python::list FindOverlappingGenes(const string& chromosome, int start, int end)
	{
		vector<string> genes;
		
		if (mGeneIntervalTrees.find(chromosome) != mGeneIntervalTrees.end())
		{
			mGeneIntervalTrees[chromosome].FindOverlapping(start, end, genes);
		}
		
		return convert_vector(genes);
	}
	
	python::list FindContainedGenes(const string& chromosome, int start, int end)
	{
		vector<string> genes;
		
		if (mGeneIntervalTrees.find(chromosome) != mGeneIntervalTrees.end())
		{
			mGeneIntervalTrees[chromosome].FindContained(start, end, genes);
		}
		
		return convert_vector(genes);
	}
	
	python::list FindNearestGenes(const string& chromosome, int position)
	{
		vector<string> genes;
		
		if (mGeneIntervalTrees.find(chromosome) != mGeneIntervalTrees.end())
		{
			mGeneIntervalTrees[chromosome].FindNearest(position, genes);
		}
		
		return convert_vector(genes);
	}
	
	string CalculateGeneLocation(const string& geneID, int position)
	{
		const Gene& gene = mGenes[geneID];
		
		if ((position < gene.start && gene.strand == "+") || (position > gene.end && gene.strand == "-"))
		{
			return "upstream";
		}
		
		if ((position > gene.end && gene.strand == "+") || (position < gene.start && gene.strand == "-"))
		{
			return "downstream";
		}
		
		bool exon = false;
		bool cds = false;
		bool utr5p = false;
		bool utr3p = false;
		for (unordered_set<string>::const_iterator transcriptIter = mGeneTranscripts[geneID].begin(); transcriptIter != mGeneTranscripts[geneID].end(); transcriptIter++)
		{
			const string& transcriptID = *transcriptIter;
			
			for (vector<Region>::const_iterator exonIter = mExons[transcriptID].begin(); exonIter != mExons[transcriptID].end(); exonIter++)
			{
				if (position >= exonIter->start && position <= exonIter->end)
				{
					exon = true;
				}
			}
			
			for (vector<Region>::const_iterator cdsIter = mCDSs[transcriptID].begin(); cdsIter != mCDSs[transcriptID].end(); cdsIter++)
			{
				if (position >= cdsIter->start && position <= cdsIter->end)
				{
					cds = true;
				}
			}
			
			if (exon && !cds && mStartCodon.find(transcriptID) != mStartCodon.end() && mStopCodon.find(transcriptID) != mStopCodon.end())
			{
				if (gene.strand == "+")
				{
					if (position < mStartCodon[transcriptID].start)
					{
						utr5p = true;
					}
					else if (position > mStopCodon[transcriptID].end)
					{
						utr3p = true;
					}
				}
				else
				{
					if (position > mStartCodon[transcriptID].end)
					{
						utr5p = true;
					}
					else if (position < mStopCodon[transcriptID].start)
					{
						utr3p = true;
					}
				}
			}
		}
		
		if (cds)
		{
			return "coding";
		}
		else if (utr5p)
		{
			return "utr5p";
		}
		else if (utr3p)
		{
			return "utr3p";
		}
		else if (exon)
		{
			return "utr";
		}
		else
		{
			return "intron";
		}
	}
	
	int CalculateGenomicPosition(const string& transcriptID, int position)
	{
		const string& geneID = mTranscriptGene[transcriptID];
		const Gene& gene = mGenes[geneID];
		const vector<Region>& exons = mExons[transcriptID];
		
		if (gene.strand == "-")
		{
			position = mTranscriptLength[transcriptID] - position + 1;
		}
		
		if (position < 1)
		{
			return exons.front().start + position - 1;
		}
		
		int localOffset = 0;
		for (vector<Region>::const_iterator exonIter = exons.begin(); exonIter != exons.end(); exonIter++)
		{
			int exonSize = exonIter->end - exonIter->start + 1;
			
			if (position <= localOffset + exonSize)
			{
				return position - localOffset - 1 + exonIter->start;
			}
			
			localOffset += exonSize;
		}
		
		return position - localOffset + exons.back().end;
	}
	
	python::list CalculateGenomicRegions(const string& transcriptID, int start, int end)
	{
		const string& geneID = mTranscriptGene[transcriptID];
		const Gene& gene = mGenes[geneID];
		const vector<Region>& exons = mExons[transcriptID];
		int transcriptLength = mTranscriptLength[transcriptID];
		
		vector<Region> regions;
		
		if (gene.strand == "-")
		{
			int tmpStart = transcriptLength - end + 1;
			int tmpEnd = transcriptLength - start + 1;
			
			start = tmpStart;
			end = tmpEnd;
		}
		
		if (start < 1)
		{
			regions.push_back(Region(exons.front().start, exons.front().start));
			return convert_vector(regions);
		}
		
		if (end > transcriptLength)
		{
			regions.push_back(Region(exons.back().end, exons.back().end));
			return convert_vector(regions);
		}
		
		int localOffset = 0;
		for (vector<Region>::const_iterator exonIter = exons.begin(); exonIter != exons.end(); exonIter++)
		{
			int localStart = start - localOffset;
			int localEnd = end - localOffset;
			
			int overlapStart = max(1, localStart) + exonIter->start - 1;
			int overlapEnd = min(exonIter->GetLength(), localEnd) + exonIter->start - 1;
			
			if (overlapStart <= overlapEnd)
			{
				regions.push_back(Region(overlapStart, overlapEnd));
			}
			
			localOffset += exonIter->GetLength();
		}
		
		return convert_vector(regions);
	}
	
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & mGenes;
		ar & mGeneTranscripts;
		ar & mTranscriptGene;
		ar & mTranscriptLength;
		ar & mExons;
		ar & mCDSs;
		ar & mStartCodon;
		ar & mStopCodon;
		ar & mGeneIntervalTrees;
	}
	
private:
	unordered_map<string,Gene> mGenes;
	unordered_map<string,unordered_set<string> > mGeneTranscripts;
	unordered_map<string,string> mTranscriptGene;
	unordered_map<string,int> mTranscriptLength;
	unordered_map<string,vector<Region> > mExons;
	unordered_map<string,vector<Region> > mCDSs;
	unordered_map<string,Region> mStartCodon;
	unordered_map<string,Region> mStopCodon;
	unordered_map<string,IntervalTree<string> > mGeneIntervalTrees;
};

BOOST_CLASS_VERSION(GeneModels, 1)

class PyIntervalTree
{
public:
	PyIntervalTree(const python::list& pyIntervals)
	{
		vector<Interval<int> > intervals;
		
		for (int i = 0; i < len(pyIntervals); ++i)
		{
			const python::tuple& pyInterval = python::extract<python::tuple>(pyIntervals[i]);
			if (len(pyInterval) != 3)
			{
				throw std::invalid_argument("Interval requires id, start, stop");
			}
			
			Interval<int> interval;
			interval.value = python::extract<int>(pyInterval[0]);
			interval.start = python::extract<int>(pyInterval[1]);
			interval.stop = python::extract<int>(pyInterval[2]);
			
			intervals.push_back(interval);
		}
		
		mIntervalTree = IntervalTree<int>(intervals);
	}
	
	python::list FindOverlapping(int start, int stop)
	{
		vector<int> overlapping;
		mIntervalTree.FindOverlapping(start, stop, overlapping);
		return convert_vector(overlapping);
	}
	
	python::list FindContained(int start, int stop)
	{
		vector<int> contained;
		mIntervalTree.FindContained(start, stop, contained);
		return convert_vector(contained);
	}
	
	python::list FindNearest(int position)
	{
		vector<int> nearest;
		mIntervalTree.FindNearest(position, nearest);
		return convert_vector(nearest);
	}
	
private:
	IntervalTree<int> mIntervalTree;
};


BOOST_PYTHON_MODULE(pygenes)
{
	using namespace python;
	
	class_<GeneModels>("GeneModels")
		.def("load_ensembl_gtf", &GeneModels::LoadEnsemblGTF)
		.def("load_binary", &GeneModels::LoadBinary)
		.def("save_binary", &GeneModels::SaveBinary)
		.def("get_gene", &GeneModels::GetGene)
		.def("get_transcript_gene", &GeneModels::GetTranscriptGene)
		.def("find_overlapping_genes", &GeneModels::FindOverlappingGenes)
		.def("find_contained_genes", &GeneModels::FindContainedGenes)
		.def("find_nearest_genes", &GeneModels::FindNearestGenes)
		.def("calculate_gene_location", &GeneModels::CalculateGeneLocation)
		.def("calculate_genomic_position", &GeneModels::CalculateGenomicPosition)
		.def("calculate_genomic_regions", &GeneModels::CalculateGenomicRegions)
	;
	
	class_<Region>("Region")
		.add_property("start", &Region::GetStart)
		.add_property("end", &Region::GetEnd)
		.add_property("length", &Region::GetLength)
	;
	
	class_<Gene>("Gene")
		.add_property("id", &Gene::GetID)
		.add_property("name", &Gene::GetName)
		.add_property("source", &Gene::GetSource)
		.add_property("chromosome", &Gene::GetChromosome)
		.add_property("strand", &Gene::GetStrand)
		.add_property("start", &Gene::GetStart)
		.add_property("end", &Gene::GetEnd)
		.add_property("length", &Gene::GetLength)
	;
	
	class_<PyIntervalTree>("IntervalTree", init<python::list>())
		.def("find_overlapping", &PyIntervalTree::FindOverlapping)
		.def("find_contained", &PyIntervalTree::FindContained)
		.def("find_nearest", &PyIntervalTree::FindNearest)
	;
}

