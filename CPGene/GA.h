/** --------------------------------------------------------------------------------------------
	 G³ówna klasa algorytmu genetycznego GA 
	 GA jest rozszerzana przez klasy pochodne specyficzne dla 
	 swoich typów obs³ugiwanych chromosomów GAFloat, GAString. 
	 Zawiera metody wspólne dla wszystkich klas pochodnych.
	 ---------------------------------------------------------------------------------------------*/
#include <vector>
#include "Chromosome.h"
#pragma once
class GA
{
public:
	/** rozmiar chromosomu (iloœæ genów) */
	int chromosomeDim;

	/** liczba chomosomów w populacji */
	int populationDim;

	/** rodzaj krzy¿owania */
	int crossoverType = 0;

	/** prawdopodobieñstwo wyst¹pienia mutacji podczas parowania chromosomów  */
	double mutationProb;

	/** limit pokoleñ */
	int maxGenerations;

	/** liczba wstêpnych uruchomieñ */
	int numPrelimRuns;

	/** max pokoleñ we wstêpnych uruchomieniach, pozwala na stworzenie lepszej puli chromosomów startowych */
	int maxPrelimGenerations;

	/** procentowa szansa mutacji losowej wystêpuj¹cej niezale¿nie od rankingu dopasowania chromosomów */
	int randomSelectionChance;

	/** prawdopodobieñstwo krzy¿owania podczas parowania */
	double crossoverProb;

	/** indeks najlepiej dopasowanego (fit) chromosomu populacji */
	int bestFitnessChromIndex = 0;

	/** indeks najgorzej dopasowanego (fit) chromosomu populacji */
	int worstFitnessChromIndex = 0;

	/** wektor chromosomów aktualnego pokolenia */
	std::vector<Chromosome*> chromosomes;

	/** tymczasowey wektor chromosomów nastêpnego pokolenia */
	std::vector<Chromosome*> chromNextGen;

	/** wektor puli chromosomów wstêpnych */
	std::vector<Chromosome*> prelimChrom;
	
	/** œrednie odchylenie pokolenia */
	std::vector<double> genAvgDeviation;

	/** œrednie dopasowanie pokolenia */
	std::vector<double> genAvgFitness;

	/** obliczanie statystyk ka¿dego pokolenia podczas ewolucji */
	bool computeStatistics = false;


	/** Inicjalizacja GA parametrami	 */
	GA(int chromosomeDim,
		int populationDim,
		double crossoverProb,
		int randomSelectionChance,
		int maxGenerations,
		int numPrelimRuns,
		int maxPrelimGenerations,
		double mutationProb,
		int crossoverType,
		bool computeStatistics);

	/** inicjalizacja chromosomów w sposób losowy */
	virtual void initPopulation() = 0;

	/** wykonanie losowej mutacji na wskazanym chromosomie */
	virtual void doRandomMutation(int iChromIndex) = 0;

	/** krzy¿owanie jednopunktowe na obydwu wskazanych chromosomach */
	virtual void doOnePtCrossover(Chromosome* Chrom1, Chromosome* Chrom2) = 0;

	/** kry¿owanie dwupunktowe na obydwu wskazanych chromosomach */
	virtual void doTwoPtCrossover(Chromosome* Chrom1, Chromosome* Chrom2) = 0;

	/** krzy¿owanie jednolite  */
	virtual void doUniformCrossover(Chromosome* Chrom1, Chromosome* Chrom2) = 0;

	/** zwraca wartoœæ dopasowania dla danego chromosomu */
	virtual double getFitness(int iChromIndex) = 0;

	/** Œrednie odchylenie liczone na podstawie wartoœci dopasowania wszystkich chromosomów w populacji */
	virtual double getAvgDeviationAmongChroms();
	
	/** zwraca œrednie odchylenie populacji dla wybranego pokolenia */
	virtual double getAvgDeviation(int iGeneration);

	/** zwraca œrednie dopasowanie populacji dla wybranego pokolenia */
	virtual double getAvgFitness(int iGeneration);

	/** zwraca wartoœæ prawdopodobieñstwa mutacji */
	virtual double getMutationProb();

	/** zwraca limit liczby pokoleñ	*/
	virtual int getMaxGenerations();

	/** zwraca liczbê uruchomieñ wstêpnych przed wykonaniem ewolucji */
	virtual int getNumPrelimRuns();

	/** zwraca liczbê pokoleñ dla uruchomieñ wstêpnych */
	virtual int getMaxPrelimGenerations();

	/** zwraca prawdopodobieñstwo doboru losowego */
	virtual int getRandomSelectionChance();

	/** zwraca prawdopodobieñstwo krzy¿owania */
	virtual double getCrossoverProb();

	/** zwraca rozmiar chromosomu (iloœæ genów)	*/
	virtual int getChromosomeDim();

	/** zwraca rozmiar populacji (iloœæ chromosomów) */
	virtual int getPopulationDim();

	/** zwraca wybrany typ krzy¿owania */
	virtual int getCrossoverType();

	/** zwraca informacjê o obliczaniu statystyk */
	virtual bool getComputeStatistics();

	/** zwraca najlepiej dopasowany chromosom w populacji */
	virtual Chromosome *getFittestChromosome();

	/** zwraca wartoœæ dopasowania najlepiej dopasowanego chromosomu w populacji */
	virtual double getFittestChromosomesFitness();

	/** zwraca losow¹ liczbê ca³kowit¹ z zakresu [0:parametr] */
	virtual int getRandom(int upperBound);

	/** zwraca losow¹ liczbê rzeczywist¹ z zakresu [0.0:parametr]	*/
	virtual double getRandom(double upperBound);

	/* funkcja przeprowadzaj¹ca algorytm ewolucji */
	virtual int evolve();

	/** oblicza œrednie dopasowanie chromosomów w populacji */
	virtual double getAvgFitness();

	/** wybór par rodziców z populacji dajacy pierwszeñstwo chromosomom dobrze dopasowanym	*/
	virtual void selectTwoParents(std::vector<int> &indexParents);

	/** Oblicza pozycjê wartoœci dopasowania chromosomu w rankingu pokolenia, wy¿sza wartoœæ oznacza wy¿sz¹ pozycjê*/
	virtual int getFitnessRank(double fitness);

	/** Obliczanie rankingu dopasowania dla wszystkich chromosomów populacji */
	virtual void computeFitnessRankings();

	/** Parowanie chromosomów w celu stworzenia nowego pokolenia zgodnie z zasadami pierwszeñstwa 
	chromosomów lepiej dopasowanych oraz elityzmem (zapewnienie przetrwania 2 najlepszych chromosomów) */
	virtual void doGeneticMating();

	/** Kopiowanie chromosomów z puli nastêpnego pokolenia do puli pokolenia aktualnego, mo¿liwoœæ mutacji */ 
	virtual void copyNextGenToThisGen();

	/** Przekazanie danych o chromosomach do widoku konsoli */
	virtual void addChromosomesToLog(int iGeneration, int iNumChromosomesToDisplay);

	/** konwertowanie informacji w postaci binarnej w string do liczby ca³kowitej '1101' --> 13 */
	virtual long long binaryStrToInt(const std::string &sBinary);

	~GA();
};