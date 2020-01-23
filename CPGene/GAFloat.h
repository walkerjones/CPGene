#pragma once
#include "GA.h"
#include "ChromFloat.h"
class GAFloat :
	public GA
{
public:
	/** liczba pozycji dziesi�tnych oznaczaj�ca precyzj� liczb zmiennoprzecinkowych	gdy chromosom 
	ma rozmiar 10, a precyzja dziesi�tnych wynosi 9, liczby w chromosomach maj� format "0.123456789" */
	int decPtsPrecision = 0;

	/** ograniczenie warto�ci chromosom�w do liczb nieujemnych */
	bool positiveNumOnly = false;

	/** zwraca chromosom rzutowany na ChromFloat * */
	virtual ChromFloat *getChromosome(int index);

	/** podmiana dw�ch gen�w w chromosomie zgodnie z podanym indeksem */
	void doRandomMutation(int iChromIndex) override;


	/** tworzy losowe chromosomy na podstawie dost�pnych gen�w z puli  */
	void initPopulation() override;
	
	/** jednopunktowe krzy�owanie chromosom�w */
	void doOnePtCrossover(Chromosome *Chrom1, Chromosome *Chrom2) override;

	/** podmiana warto�ci gen�w w dw�ch chromosomach na wskazanej pozycji */
	virtual void swapGene(int geneIndex, ChromFloat *chrom1, ChromFloat *chrom2);

	/** dwupunkotwe krzy�owanie chromosom�w */
	void doTwoPtCrossover(Chromosome *Chrom1, Chromosome *Chrom2) override;

	/** jednolite krzy�owanie chromosom�w  */
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

