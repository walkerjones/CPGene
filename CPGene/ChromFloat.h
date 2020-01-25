/**
 * Klasa rozszerzaj�ca chromosome dla typu chromosom�w zapisanych w postaci liczb zmiennoprzecinkowych
 */
#pragma once
#include "Chromosome.h"
#include <vector>
#include <string>

class ChromFloat :
	public Chromosome
{
public:
	/** chromosom - wektor gen�w typu double */
	std::vector<double> genes; 

	/** Konstruktor tworzy wektor dla gen�w o rozmiarze podanym w parametrze */
	ChromFloat(int iGenesDim);

	/** Zwraca wektor gen�w */
	std::vector<double> getGenes();

	/** Zwraca chromosom w postaci string */
	std::string getGenesAsStr();

	/** Zwraca sformatowany �a�cuch typu string zawieraj�cy geny */
	std::string getGenesAsString();

	/** Kopiuje geny ze wzkazanego chromosomu do obecnego, nadpisuj�c je */
	void copyChromGenes(Chromosome *chromosome);

	/** Zwraca pojedynczy gen ze wzkazanej pozycji w chromosomie */
	double getGene(int iGene);

	/** Ustawia wkazany gen na podan� warto�� */
	void setGene(int iGene, double value);

	/** Oblicza liczb� wsp�lnych gen�w dla danego chromosomu i przekazanego argumentem */
	int getNumGenesInCommon(Chromosome* chromosome);

	~ChromFloat();
};

