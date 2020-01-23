#pragma once
#include "GA.h"
#include "ChromFloat.h"
class GAFloat :
	public GA
{
public:
	/** liczba pozycji dziesiêtnych oznaczaj¹ca precyzjê liczb zmiennoprzecinkowych	gdy chromosom 
	ma rozmiar 10, a precyzja dziesiêtnych wynosi 9, liczby w chromosomach maj¹ format "0.123456789" */
	int decPtsPrecision = 0;

	/** ograniczenie wartoœci chromosomów do liczb nieujemnych */
	bool positiveNumOnly = false;

	/** zwraca chromosom rzutowany na ChromFloat * */
	virtual ChromFloat *getChromosome(int index);

	/** podmiana dwóch genów w chromosomie zgodnie z podanym indeksem */
	void doRandomMutation(int iChromIndex) override;


	/** tworzy losowe chromosomy na podstawie dostêpnych genów z puli  */
	void initPopulation() override;
	
	/** jednopunktowe krzy¿owanie chromosomów */
	void doOnePtCrossover(Chromosome *Chrom1, Chromosome *Chrom2) override;

	/** podmiana wartoœci genów w dwóch chromosomach na wskazanej pozycji */
	virtual void swapGene(int geneIndex, ChromFloat *chrom1, ChromFloat *chrom2);

	/** dwupunkotwe krzy¿owanie chromosomów */
	void doTwoPtCrossover(Chromosome *Chrom1, Chromosome *Chrom2) override;

	/** jednolite krzy¿owanie chromosomów  */
	void doUniformCrossover(Chromosome *Chrom1, Chromosome *Chrom2) override;

	/** Inicjalizacja obiektu z parametrami */

	GAFloat(int chromosomeDim, 
		int populationDim,
		double crossoverProb, 
		int randomSelectionChance, 
		int maxGenerations, 
		int numPrelimRuns, 
		int maxPrelimGenerations, 
		double mutationProb, 
		int crossoverType, 
		int decPtsPrecision, 
		bool positiveNumOnly, 
		bool computeStatistics);


	~GAFloat();
	double getFitness(int iChromIndex);

#ifdef CurveFit
	public:	
	int curveDim = 0;
	std::vector<double> curveData;
	virtual void setCurveData(std::vector<double> &CurveData);
#endif
};

