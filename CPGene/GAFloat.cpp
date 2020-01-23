#include "PCH.h"




GAFloat::GAFloat(int chromosomeDim, int populationDim, double crossoverProb, int randomSelectionChance, int maxGenerations, int numPrelimRuns, int maxPrelimGenerations, double mutationProb, int crossoverType, int decPtsPrecision, bool positiveNumOnly, bool computeStatistics) :
	GA(chromosomeDim, populationDim, crossoverProb, randomSelectionChance, maxGenerations, numPrelimRuns, maxPrelimGenerations, mutationProb, crossoverType, computeStatistics)
{

	if (decPtsPrecision < 0)
	{
		std::cout << "Precyzja nie moze byc ujemna" << std::endl;
	}
	if (chromosomeDim < 1)
	{
		std::cout << "Chromosom musi mieæ rozmiar > 0" << std::endl;
	}

	for (int i = 0; i < populationDim; i++)
	{
		this->chromosomes[i] = new ChromFloat(chromosomeDim);
		this->chromNextGen[i] = new ChromFloat(chromosomeDim);
		this->prelimChrom[i] = new ChromFloat(chromosomeDim);
	}

	this->decPtsPrecision = decPtsPrecision;
	this->positiveNumOnly = positiveNumOnly;
 
}

GAFloat::~GAFloat()
{
}

ChromFloat * GAFloat::getChromosome(int index)
{
	return (static_cast<ChromFloat*>(this->chromosomes[index]));
}

void GAFloat::doRandomMutation(int iChromIndex)
{
	int iGene;
	double rNewGene;

	iGene = getRandom(chromosomeDim-1);

	rNewGene = (static_cast<ChromFloat*>(this->chromosomes[iChromIndex]))->genes[iGene];
	if (getRandom(100) > 50)
	{
		rNewGene = rNewGene + (rNewGene * getRandom(1000.0) / 1000.0);
	}
	else
	{
		rNewGene = rNewGene - (rNewGene * getRandom(1000.0) / 1000.0);
	}

	(static_cast<ChromFloat*>(this->chromosomes[iChromIndex]))->genes[iGene] = rNewGene;
}

void GAFloat::initPopulation()
{
	for (int i = 0; i < populationDim; i++)
	{
		for (int iGene = 0; iGene < chromosomeDim; iGene++)
		{
			if ((this->positiveNumOnly == true) || (getRandom(100) > 50))
			{
				(static_cast<ChromFloat*>(this->chromosomes[i]))->genes[iGene] = getRandom(1000.0) / (1 + getRandom(1000.0));
			}
			else
			{
				(static_cast<ChromFloat*>(this->chromosomes[i]))->genes[iGene] = -getRandom(1000.0) / (1 + getRandom(1000.0));
			}
		}
		this->chromosomes[i]->fitness = getFitness(i);
	}
}

void GAFloat::doOnePtCrossover(Chromosome * Chrom1, Chromosome * Chrom2)
{
	int crossoverPt = getRandom(chromosomeDim-1);
	swapGene(crossoverPt, static_cast<ChromFloat*>(Chrom1), static_cast<ChromFloat*>(Chrom2));
}

void GAFloat::swapGene(int geneIndex, ChromFloat * chrom1, ChromFloat * chrom2)
{
	double temp = chrom1->genes[geneIndex];
	chrom1->genes[geneIndex] = chrom2->genes[geneIndex];
	chrom2->genes[geneIndex] = temp;
}

void GAFloat::doTwoPtCrossover(Chromosome * Chrom1, Chromosome * Chrom2)
{
	int crossPt1 = getRandom(chromosomeDim-1);
	int crossPt2 = getRandom(chromosomeDim-1);

	if (crossPt1 != crossPt2)
	{
		for (int geneIndex = crossPt1; geneIndex < crossPt2; geneIndex++)
		{
			swapGene(geneIndex, static_cast<ChromFloat*>(Chrom1), static_cast<ChromFloat*>(Chrom2));
		}
	}
}

void GAFloat::doUniformCrossover(Chromosome * Chrom1, Chromosome * Chrom2)
{
	std::string sGene1, sGene2;
	std::string sNewChrom1 = "", sNewChrom2 = "";

	for (int geneIndex = 0; geneIndex < chromosomeDim; geneIndex++)
	{
		if (getRandom(100) > 50)
		{
			swapGene(geneIndex, static_cast<ChromFloat*>(Chrom1), static_cast<ChromFloat*>(Chrom2));
		}
	}
}
