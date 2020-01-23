#include "PCH.h"


GAString::~GAString()
{
}

ChromChars * GAString::getChromosome(int index)
{
	return (static_cast<ChromChars*>(this->chromosomes[index]));
}

GAString::GAString(
	int chromosomeDim, 
	int populationDim, 
	double crossoverProb, 
	int randomSelectionChance, 
	int maxGenerations, 
	int numPrelimRuns, 
	int maxPrelimGenerations, 
	double mutationProb, 
	int chromDecPts, 
	const std::string & possGeneValues, 
	int crossoverType, 
	bool computeStatistics) : 
	GA(chromosomeDim, 
		populationDim, 
		crossoverProb, 
		randomSelectionChance, 
		maxGenerations, 
		numPrelimRuns, 
		maxPrelimGenerations, 
		mutationProb, 
		crossoverType, 
		computeStatistics)
{
	if (possGeneValues.length() < 2)
	{
		std::cout << "Pula mozliwych genow musi wynosic minimum 2" << std::endl;
	}

	this->chromDecPts = chromDecPts;
	this->possGeneValues = possGeneValues;

	for (int i = 0; i < populationDim; i++)
	{
		this->chromosomes[i] = new ChromChars(chromosomeDim);
		this->chromNextGen[i] = new ChromChars(chromosomeDim);
		this->prelimChrom[i] = new ChromChars(chromosomeDim);
	}
	
}

double GAString::chromStrToFloat(const std::string & sChromosome, int iNumDecPts)
{
	std::string sFloat;
	int iLen;

	if (iNumDecPts == 0)
	{
		return static_cast<double>(binaryStrToInt(sChromosome));
	}
	else
	{
		iLen = sChromosome.length() - iNumDecPts;
		sFloat = sChromosome.substr(0, iLen) + "." + sChromosome.substr(iLen, iNumDecPts);
		return (std::stod(sFloat));
	}
}

double GAString::getChromValAsDouble(const std::string & sChromosome)
{
	return (chromStrToFloat(sChromosome, this->chromDecPts));
}

char GAString::getRandomGeneFromPossGenes()
{
	int iRandomIndex = getRandom((static_cast<int>((this->possGeneValues.length())-1)));
	return (this->possGeneValues[iRandomIndex]);
}

void GAString::doRandomMutation(int iChromIndex)
{
	int iGene1, iGene2;
	char cTemp;

	iGene1 = getRandom(chromosomeDim-1);
	iGene2 = getRandom(chromosomeDim-1);

	cTemp = (static_cast<ChromChars*>(this->chromosomes[iChromIndex]))->genes[iGene1];
	(static_cast<ChromChars*>(this->chromosomes[iChromIndex]))->genes[iGene1] = (static_cast<ChromChars*>(this->chromosomes[iChromIndex]))->genes[iGene2];
	(static_cast<ChromChars*>(this->chromosomes[iChromIndex]))->genes[iGene2] = cTemp;
}

void GAString::initPopulation()
{
	for (int i = 0; i < populationDim; i++)
	{
		for (int iGene = 0; iGene < chromosomeDim; iGene++)
		{
			(static_cast<ChromChars*>(this->chromosomes[i]))->genes[iGene] = getRandomGeneFromPossGenes();
		}
		this->chromosomes[i]->fitness = this->getFitness(i);
	}
}

void GAString::doOnePtCrossover(Chromosome * Chrom1, Chromosome * Chrom2)
{
	std::string sNewChrom1, sNewChrom2;
	int iCrossoverPoint;
	std::string sChrom1, sChrom2;

	iCrossoverPoint = getRandom(chromosomeDim - 2);
	sChrom1 = Chrom1->getGenesAsStr();
	sChrom2 = Chrom2->getGenesAsStr();

	sNewChrom1 = sChrom1.substr(0, iCrossoverPoint) + sChrom2.substr(iCrossoverPoint, chromosomeDim - iCrossoverPoint);
	sNewChrom2 = sChrom2.substr(0, iCrossoverPoint) + sChrom1.substr(iCrossoverPoint, chromosomeDim - iCrossoverPoint);

	(static_cast<ChromChars*>(Chrom1))->setGenesFromStr(sNewChrom1);
	(static_cast<ChromChars*>(Chrom2))->setGenesFromStr(sNewChrom2);
}

void GAString::doTwoPtCrossover(Chromosome * Chrom1, Chromosome * Chrom2)
{
	std::string sNewChrom1, sNewChrom2;
	int iCrossoverPoint1, iCrossoverPoint2;
	std::string sChrom1, sChrom2;

	iCrossoverPoint1 = 1 + getRandom(chromosomeDim - 2);
	iCrossoverPoint2 = iCrossoverPoint1 + 1 + getRandom(chromosomeDim - iCrossoverPoint1 - 1);

	if (iCrossoverPoint2 == (iCrossoverPoint1 + 1))
	{
		doOnePtCrossover(Chrom1, Chrom2);
	}
	else
	{
		sChrom1 = Chrom1->getGenesAsStr();
		sChrom2 = Chrom2->getGenesAsStr();
		sNewChrom1 = sChrom1.substr(0, iCrossoverPoint1) + sChrom2.substr(iCrossoverPoint1, iCrossoverPoint2 - iCrossoverPoint1) + sChrom1.substr(iCrossoverPoint2, chromosomeDim - iCrossoverPoint2);
		sNewChrom2 = sChrom2.substr(0, iCrossoverPoint1) + sChrom1.substr(iCrossoverPoint1, iCrossoverPoint2 - iCrossoverPoint1) + sChrom2.substr(iCrossoverPoint2, chromosomeDim - iCrossoverPoint2);

		(static_cast<ChromChars*>(Chrom1))->setGenesFromStr(sNewChrom1);
		(static_cast<ChromChars*>(Chrom2))->setGenesFromStr(sNewChrom2);
	}
}

void GAString::doUniformCrossover(Chromosome * Chrom1, Chromosome * Chrom2)
{
	int iGeneToSwap;
	char cGene;

	StringBuilder *sbChrom1 = new StringBuilder(Chrom1->getGenesAsStr());
	StringBuilder *sbChrom2 = new StringBuilder(Chrom2->getGenesAsStr());

	for (int i = 0; i < chromosomeDim; i++)
	{
		if (getRandom(100) > 50)
		{
			iGeneToSwap = getRandom(chromosomeDim-1);
			cGene = sbChrom1->charAt(iGeneToSwap);

			sbChrom1->setCharAt(iGeneToSwap, sbChrom2->charAt(iGeneToSwap));
			sbChrom2->setCharAt(iGeneToSwap, cGene);
		}
	}

	(static_cast<ChromChars*>(Chrom1))->setGenesFromStr(sbChrom1->toString());
	(static_cast<ChromChars*>(Chrom2))->setGenesFromStr(sbChrom2->toString());

	delete sbChrom2;
	delete sbChrom1;
}