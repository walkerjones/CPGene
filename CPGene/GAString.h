#pragma once

class GAString :
	public GA
{
public:
	/** precyzja dziesiêtna, liczba pól po przecinku jeœli chromosom przedstawia liczbê rzeczywist¹ w postaci string */
	int chromDecPts = 0;

	/** pula mo¿liwych postaci genów  */
	std::string possGeneValues;

	/** zwraca chromosom rzutowany na ChromChars * */
	ChromChars *getChromosome(int index);

	/** Utworzenie obiektu z parametrami */
	GAString(int chromosomeDim,
		int populationDim,
		double crossoverProb,
		int randomSelectionChance,
		int maxGenerations,
		int numPrelimRuns,
		int maxPrelimGenerations,
		double mutationProb,
		int chromDecPts,
		const std::string &possGeneValues,
		int crossoverType,
		bool computeStatistics);

	/** Konwersja chromosomu w postaci string do liczby rzeczywistej z uwzglêdnieniem liczby miejsc po przecinku */
	double chromStrToFloat(const std::string &sChromosome, int iNumDecPts);

	/** KOnwersja chromosomu w postaci string do liczby rzeczywistej */
	double getChromValAsDouble(const std::string &sChromosome);

	/** Zwraca losow¹ wartoœæ genu z puli mo¿liwych wartoœci */
	char getRandomGeneFromPossGenes();

	/** Losowo podmienia wartoœci genów w dwóch chromosomach na podanej pozycji */
	void doRandomMutation(int iChromIndex);

	/** tworzy losowe chromosomy */
	void initPopulation();

	/** jednopunktowe krzy¿owanie */
	void doOnePtCrossover(Chromosome *Chrom1, Chromosome *Chrom2);

	/** dwupunktowe krzy¿owanie */
	void doTwoPtCrossover(Chromosome *Chrom1, Chromosome *Chrom2);

	/** jednolite krzy¿owanie */
	void doUniformCrossover(Chromosome *Chrom1, Chromosome *Chrom2);

	~GAString();
	virtual double getFitness(int iChromIndex);

#ifdef Knapsack
	double weight, value, limit;
	std::vector<double> knapWeight;
	std::vector<double> knapValue;
	void printResult(std::ofstream &txtout);
#endif
};