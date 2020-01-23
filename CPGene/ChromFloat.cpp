#include "PCH.h"


ChromFloat::ChromFloat(int iGenesDim)
{
	genes = std::vector<double>(iGenesDim);
}

std::vector<double> ChromFloat::getGenes()
{
	return genes;
}

std::string ChromFloat::toString()
{
	return getGenesAsString();
}

int ChromFloat::getNumGenesInCommon(Chromosome * chromosome)
{
	return genes.size();
}

std::string ChromFloat::getGenesAsString()
{
	std::string sGenes = "";
	for (size_t i = 0; i < genes.size(); i++)
	{
		sGenes += " " + std::to_string(genes[i]) + ",";
	}
	return sGenes;
}

void ChromFloat::copyChromGenes(Chromosome * chromosome)
{
	ChromFloat *chromFloat = static_cast<ChromFloat*>(chromosome);
	for (size_t iGene = 0; iGene < genes.size(); iGene++)
	{
		this->genes[iGene] = chromFloat->genes[iGene];
	}
}

double ChromFloat::getGene(int iGene)
{
	return genes[iGene];
}

void ChromFloat::setGene(int iGene, double value)
{
	genes[iGene] = value;
}

std::string ChromFloat::getGenesAsStr()
{
	std::string sGenes = "";
	for (size_t i = 0; i < genes.size(); i++)
	{
		sGenes += std::to_string(genes[i]) + ",";
	}
	return sGenes;
}

ChromFloat::~ChromFloat()
{
}
