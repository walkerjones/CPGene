#include "PCH.h"


ChromChars::ChromChars(int iGenesDim)
{
	genes = std::vector<char>(iGenesDim);
}

int ChromChars::getNumGenesInCommon(Chromosome * chromosome)
{
	int numGenesInCommon = 0;
	std::string chromGenes = chromosome->getGenesAsStr();

	for (size_t i = 0; i < genes.size(); i++)
	{
		if (this->genes[i] == chromGenes[i])
		{
			numGenesInCommon++;
		}
	}
	return (numGenesInCommon);
}

std::vector<char> ChromChars::getGenes()
{
	return (genes);
}

std::string ChromChars::getGenesAsStr()
{
	std::string sGenes = "";

	for (size_t i = 0; i < genes.size(); i++)
	{
		sGenes += genes[i];
	}
	return (sGenes);
}

char ChromChars::getGene(int iGene)
{
	return (this->genes[iGene]);
}

void ChromChars::setGenesFromStr(const std::string & sChromosome)
{
	for (size_t i = 0; i < genes.size(); i++)
	{
		this->genes[i] = sChromosome[i];
	}
}


void ChromChars::copyChromGenes(Chromosome * chromosome)
{
	ChromChars *chromChars = static_cast<ChromChars*>(chromosome);
	for (size_t iGene = 0; iGene < genes.size(); iGene++)
	{
		this->genes[iGene] = chromChars->genes[iGene];
	}
}

ChromChars::~ChromChars()
{
}
