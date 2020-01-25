/** --------------------------------------------------------------------------------------------
	 G��wna klasa algorytmu genetycznego GA 
	 GA jest rozszerzana przez klasy pochodne specyficzne dla 
	 swoich typ�w obs�ugiwanych chromosom�w GAFloat, GAString. 
	 Zawiera metody wsp�lne dla wszystkich klas pochodnych.
	 ---------------------------------------------------------------------------------------------*/
#include <vector>
#include "Chromosome.h"
#pragma once
class GA
{
public:
	/** rozmiar chromosomu (ilo�� gen�w) */
	int chromosomeDim;

	/** liczba chomosom�w w populacji */
	int populationDim;

	/** rodzaj krzy�owania */
	int crossoverType = 0;

	/** prawdopodobie�stwo wyst�pienia mutacji podczas parowania chromosom�w  */
	double mutationProb;

	/** limit pokole� */
	int maxGenerations;

	/** liczba wst�pnych uruchomie� */
	int numPrelimRuns;

	/** max pokole� we wst�pnych uruchomieniach, pozwala na stworzenie lepszej puli chromosom�w startowych */
	int maxPrelimGenerations;

	/** procentowa szansa mutacji losowej wyst�puj�cej niezale�nie od rankingu dopasowania chromosom�w */
	int randomSelectionChance;

	/** prawdopodobie�stwo krzy�owania podczas parowania */
	double crossoverProb;

	/** indeks najlepiej dopasowanego (fit) chromosomu populacji */
	int bestFitnessChromIndex = 0;

	/** indeks najgorzej dopasowanego (fit) chromosomu populacji */
	int worstFitnessChromIndex = 0;

	/** wektor chromosom�w aktualnego pokolenia */
	std::vector<Chromosome*> chromosomes;

	/** tymczasowey wektor chromosom�w nast�pnego pokolenia */
	std::vector<Chromosome*> chromNextGen;

	/** wektor puli chromosom�w wst�pnych */
	std::vector<Chromosome*> prelimChrom;
	
	/** �rednie odchylenie pokolenia */
	std::vector<double> genAvgDeviation;

	/** �rednie dopasowanie pokolenia */
	std::vector<double> genAvgFitness;

	/** obliczanie statystyk ka�dego pokolenia podczas ewolucji */
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

	/** inicjalizacja chromosom�w w spos�b losowy */
	virtual void initPopulation() = 0;

	/** wykonanie losowej mutacji na wskazanym chromosomie */
	virtual void doRandomMutation(int iChromIndex) = 0;

	/** krzy�owanie jednopunktowe na obydwu wskazanych chromosomach */
	virtual void doOnePtCrossover(Chromosome* Chrom1, Chromosome* Chrom2) = 0;

	/** kry�owanie dwupunktowe na obydwu wskazanych chromosomach */
	virtual void doTwoPtCrossover(Chromosome* Chrom1, Chromosome* Chrom2) = 0;

	/** krzy�owanie jednolite  */
	virtual void doUniformCrossover(Chromosome* Chrom1, Chromosome* Chrom2) = 0;

	/** zwraca warto�� dopasowania dla danego chromosomu */
	virtual double getFitness(int iChromIndex) = 0;

	/** �rednie odchylenie liczone na podstawie warto�ci dopasowania wszystkich chromosom�w w populacji */
	virtual double getAvgDeviationAmongChroms();
	
	/** zwraca �rednie odchylenie populacji dla wybranego pokolenia */
	virtual double getAvgDeviation(int iGeneration);

	/** zwraca �rednie dopasowanie populacji dla wybranego pokolenia */
	virtual double getAvgFitness(int iGeneration);

	/** zwraca warto�� prawdopodobie�stwa mutacji */
	virtual double getMutationProb();

	/** zwraca limit liczby pokole�	*/
	virtual int getMaxGenerations();

	/** zwraca liczb� uruchomie� wst�pnych przed wykonaniem ewolucji */
	virtual int getNumPrelimRuns();

	/** zwraca liczb� pokole� dla uruchomie� wst�pnych */
	virtual int getMaxPrelimGenerations();

	/** zwraca prawdopodobie�stwo doboru losowego */
	virtual int getRandomSelectionChance();

	/** zwraca prawdopodobie�stwo krzy�owania */
	virtual double getCrossoverProb();

	/** zwraca rozmiar chromosomu (ilo�� gen�w)	*/
	virtual int getChromosomeDim();

	/** zwraca rozmiar populacji (ilo�� chromosom�w) */
	virtual int getPopulationDim();

	/** zwraca wybrany typ krzy�owania */
	virtual int getCrossoverType();

	/** zwraca informacj� o obliczaniu statystyk */
	virtual bool getComputeStatistics();

	/** zwraca najlepiej dopasowany chromosom w populacji */
	virtual Chromosome *getFittestChromosome();

	/** zwraca warto�� dopasowania najlepiej dopasowanego chromosomu w populacji */
	virtual double getFittestChromosomesFitness();

	/** zwraca losow� liczb� ca�kowit� z zakresu [0:parametr] */
	virtual int getRandom(int upperBound);

	/** zwraca losow� liczb� rzeczywist� z zakresu [0.0:parametr]	*/
	virtual double getRandom(double upperBound);

	/* funkcja przeprowadzaj�ca algorytm ewolucji */
	virtual int evolve();

	/** oblicza �rednie dopasowanie chromosom�w w populacji */
	virtual double getAvgFitness();

	/** wyb�r par rodzic�w z populacji dajacy pierwsze�stwo chromosomom dobrze dopasowanym	*/
	virtual void selectTwoParents(std::vector<int> &indexParents);

	/** Oblicza pozycj� warto�ci dopasowania chromosomu w rankingu pokolenia, wy�sza warto�� oznacza wy�sz� pozycj�*/
	virtual int getFitnessRank(double fitness);

	/** Obliczanie rankingu dopasowania dla wszystkich chromosom�w populacji */
	virtual void computeFitnessRankings();

	/** Parowanie chromosom�w w celu stworzenia nowego pokolenia zgodnie z zasadami pierwsze�stwa 
	chromosom�w lepiej dopasowanych oraz elityzmem (zapewnienie przetrwania 2 najlepszych chromosom�w) */
	virtual void doGeneticMating();

	/** Kopiowanie chromosom�w z puli nast�pnego pokolenia do puli pokolenia aktualnego, mo�liwo�� mutacji */ 
	virtual void copyNextGenToThisGen();

	/** Przekazanie danych o chromosomach do widoku konsoli */
	virtual void addChromosomesToLog(int iGeneration, int iNumChromosomesToDisplay);

	/** konwertowanie informacji w postaci binarnej w string do liczby ca�kowitej '1101' --> 13 */
	virtual long long binaryStrToInt(const std::string &sBinary);

	~GA();
};