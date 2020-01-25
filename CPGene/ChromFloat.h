/**
 * Klasa rozszerzaj¹ca chromosome dla typu chromosomów zapisanych w postaci liczb zmiennoprzecinkowych
 */
#pragma once
#include "Chromosome.h"
#include <vector>
#include <string>

class ChromFloat :
	public Chromosome
{
public:
	/** chromosom - wektor genów typu double */
	std::vector<double> genes; 

	/** Konstruktor tworzy wektor dla genów o rozmiarze podanym w parametrze */
	ChromFloat(int iGenesDim);

	/** Zwraca wektor genów */
	std::vector<double> getGenes();

	/** Zwraca chromosom w postaci string */
	std::string getGenesAsStr();

	/** Zwraca sformatowany ³añcuch typu string zawieraj¹cy geny */
	std::string getGenesAsString();

	/** Kopiuje geny ze wzkazanego chromosomu do obecnego, nadpisuj¹c je */
	void copyChromGenes(Chromosome *chromosome);

	/** Zwraca pojedynczy gen ze wzkazanej pozycji w chromosomie */
	double getGene(int iGene);

	/** Ustawia wkazany gen na podan¹ wartoœæ */
	void setGene(int iGene, double value);

	/** Oblicza liczbê wspólnych genów dla danego chromosomu i przekazanego argumentem */
	int getNumGenesInCommon(Chromosome* chromosome);

	~ChromFloat();
};

