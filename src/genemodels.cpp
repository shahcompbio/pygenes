
#include "IntervalTree/IntervalTree.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>

using namespace std;


template<typename T>
void Print(const vector<T>& s)
{
	for (typename vector<T>::const_iterator Iter = s.begin(); Iter != s.end(); Iter++)
	{
		cout << *Iter << endl;
	}	
}


struct CRegion
{
	CRegion() : start(-1), end(-1) {}
	CRegion(int start, int end) : start(start), end(end) {}
	
	int start;
	int end;
	
	int GetStart() const { return start; }
	int GetEnd() const { return end; }
	int GetLength() const { return end - start + 1; }
	
	bool operator==(const CRegion& other) const
	{
		return start == other.start && end == other.end;
	}
	
	bool operator<(const CRegion& other) const
	{
		return start < other.start;
	}
};


struct CGene : public CRegion
{
	string id;
	string name;
	string source;
	string chromosome;
	string strand;
	
	bool operator==(const CGene& other) const
	{
		return id == other.id;
	}
};


void splitbytoken(vector<string>& fields, const string& line, char token) 
{
	stringstream linestream(line);
	string intermediate;
	while(getline(linestream, intermediate, token))
	{
		fields.push_back(intermediate);
	}
}


template<typename T>
T lexical_cast(const string& field)
{
	stringstream fieldstream(field);

	T value;
	fieldstream >> value;

	if (fieldstream.good())
	{
		stringstream errorStr;
		errorStr << "Field " << field << " could not be cast";
		throw std::invalid_argument(errorStr.str());
	}

	return value;
}


string trim(string& str)
{
	string::iterator iter1;
	for (iter1 = str.begin(); iter1 != str.end(); iter1++)
	{
		if (*iter1 != ' ')
		{
			break;
		}
	}

	string::iterator iter2;
	for (iter2 = str.end(); iter2 != iter1; iter2--)
	{
		if (*iter2 != ' ')
		{
			break;
		}
	}

	return string(iter1, iter2);
}


class CGeneModels
{
public:
	void LoadEnsemblGTF(const string& gtfFilename)
	{
		map<string, vector<int> > geneRegionStart;
		map<string, vector<int> > geneRegionEnd;
		
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

			splitbytoken(gtfFields, line, '\t');

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
			splitbytoken(featureFields, gtfFields[8], ';');
			
			string geneID;
			string transcriptID;
			string geneName;
			int exonNumber = -1;
			for (vector<string>::iterator featureIter = featureFields.begin(); featureIter != featureFields.end(); featureIter++)
			{
				string feature = *featureIter;

				feature = trim(feature);

				if (feature.empty())
				{
					continue;
				}
				
				int firstSpace = feature.find_first_of(' ');
				
				string key = feature.substr(0, firstSpace);
				string value = feature.substr(firstSpace + 2, feature.size() - firstSpace - 3);
				
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
			geneRegionStart[geneID].push_back(start);
			geneRegionEnd[geneID].push_back(end);
			
			mGeneTranscripts[geneID].insert(transcriptID);
			mTranscriptGene[transcriptID] = geneID;
			
			if (featureType == "exon")
			{
				mExons[transcriptID].push_back(CRegion(start, end));
			}
			else if (featureType == "CDS")
			{
				mCDSs[transcriptID].push_back(CRegion(start, end));
			}
			else if (featureType == "start_codon")
			{
				mStartCodon[transcriptID] = CRegion(start, end);
			}
			else if (featureType == "stop_codon")
			{
				mStopCodon[transcriptID] = CRegion(start, end);
			}
		}
		
		map<string,vector<CInterval<string> > > geneIntervals;
		
		for (map<string,CGene>::const_iterator geneIter = mGenes.begin(); geneIter != mGenes.end(); geneIter++)
		{
			const string geneID = geneIter->first;
			
			int geneStart = *min_element(geneRegionStart[geneID].begin(), geneRegionStart[geneID].end());
			int geneEnd = *max_element(geneRegionEnd[geneID].begin(), geneRegionEnd[geneID].end());
			
			mGenes[geneID].start = geneStart;
			mGenes[geneID].end = geneEnd;
			
			geneIntervals[mGenes[geneID].chromosome].push_back(CInterval<string>(geneStart, geneEnd, geneID));
		}
		
		mGeneIntervalTrees.clear();
		for (map<string,vector<CInterval<string> > >::iterator chromosomeIter = geneIntervals.begin(); chromosomeIter != geneIntervals.end(); chromosomeIter++)
		{
			mGeneIntervalTrees[chromosomeIter->first] = CIntervalTree<string>(chromosomeIter->second);
		}
		
		for (map<string,vector<CRegion> >::iterator exonsIter = mExons.begin(); exonsIter != mExons.end(); exonsIter++)
		{
			sort(exonsIter->second.begin(), exonsIter->second.end());
			
			int length = 0;
			for (vector<CRegion>::const_iterator exonIter = exonsIter->second.begin(); exonIter != exonsIter->second.end(); exonIter++)
			{
				length += exonIter->end - exonIter->start + 1;
			}
			mTranscriptLength[exonsIter->first] = length;
		}
		
		for (map<string,vector<CRegion> >::iterator cdssIter = mCDSs.begin(); cdssIter != mCDSs.end(); cdssIter++)
		{
			sort(cdssIter->second.begin(), cdssIter->second.end());
		}
	}
	
	CGene GetGene(const string& geneID)
	{
		return mGenes[geneID];
	}
	
	string GetTranscriptGene(const string& transcriptID)
	{
		return mTranscriptGene[transcriptID];
	}
	
	void FindOverlappingGenes(const string& chromosome, int start, int end, vector<string>& genes)
	{
		if (mGeneIntervalTrees.find(chromosome) != mGeneIntervalTrees.end())
		{
			mGeneIntervalTrees[chromosome].FindOverlapping(start, end, genes);
		}
	}
	
	void FindContainedGenes(const string& chromosome, int start, int end, vector<string>& genes)
	{
		if (mGeneIntervalTrees.find(chromosome) != mGeneIntervalTrees.end())
		{
			mGeneIntervalTrees[chromosome].FindContained(start, end, genes);
		}
	}
	
	void FindNearestGenes(const string& chromosome, int position, vector<string>& genes)
	{
		if (mGeneIntervalTrees.find(chromosome) != mGeneIntervalTrees.end())
		{
			mGeneIntervalTrees[chromosome].FindNearest(position, genes);
		}
	}
	
	string CalculateGeneLocation(const string& geneID, int position)
	{
		const CGene& gene = mGenes[geneID];
		
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
		for (set<string>::const_iterator transcriptIter = mGeneTranscripts[geneID].begin(); transcriptIter != mGeneTranscripts[geneID].end(); transcriptIter++)
		{
			const string& transcriptID = *transcriptIter;
			
			for (vector<CRegion>::const_iterator exonIter = mExons[transcriptID].begin(); exonIter != mExons[transcriptID].end(); exonIter++)
			{
				if (position >= exonIter->start && position <= exonIter->end)
				{
					exon = true;
				}
			}
			
			for (vector<CRegion>::const_iterator cdsIter = mCDSs[transcriptID].begin(); cdsIter != mCDSs[transcriptID].end(); cdsIter++)
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
		const CGene& gene = mGenes[geneID];
		const vector<CRegion>& exons = mExons[transcriptID];
		
		if (gene.strand == "-")
		{
			position = mTranscriptLength[transcriptID] - position + 1;
		}
		
		if (position < 1)
		{
			return exons.front().start + position - 1;
		}
		
		int localOffset = 0;
		for (vector<CRegion>::const_iterator exonIter = exons.begin(); exonIter != exons.end(); exonIter++)
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
	
	void CalculateGenomicRegions(const string& transcriptID, int start, int end, vector<CRegion>& regions)
	{
		const string& geneID = mTranscriptGene[transcriptID];
		const CGene& gene = mGenes[geneID];
		const vector<CRegion>& exons = mExons[transcriptID];
		int transcriptLength = mTranscriptLength[transcriptID];
		
		if (gene.strand == "-")
		{
			int tmpStart = transcriptLength - end + 1;
			int tmpEnd = transcriptLength - start + 1;
			
			start = tmpStart;
			end = tmpEnd;
		}
		
		if (start < 1)
		{
			regions.push_back(CRegion(exons.front().start, exons.front().start));
			return;
		}
		
		if (end > transcriptLength)
		{
			regions.push_back(CRegion(exons.back().end, exons.back().end));
			return;
		}
		
		int localOffset = 0;
		for (vector<CRegion>::const_iterator exonIter = exons.begin(); exonIter != exons.end(); exonIter++)
		{
			int localStart = start - localOffset;
			int localEnd = end - localOffset;
			
			int overlapStart = max(1, localStart) + exonIter->start - 1;
			int overlapEnd = min(exonIter->GetLength(), localEnd) + exonIter->start - 1;
			
			if (overlapStart <= overlapEnd)
			{
				regions.push_back(CRegion(overlapStart, overlapEnd));
			}
			
			localOffset += exonIter->GetLength();
		}
	}
	
private:
	map<string,CGene> mGenes;
	map<string,set<string> > mGeneTranscripts;
	map<string,string> mTranscriptGene;
	map<string,int> mTranscriptLength;
	map<string,vector<CRegion> > mExons;
	map<string,vector<CRegion> > mCDSs;
	map<string,CRegion> mStartCodon;
	map<string,CRegion> mStopCodon;
	map<string,CIntervalTree<string> > mGeneIntervalTrees;
};


