/**
* Klasa rozszerza Chromosome dla chromosom�w typu Char
*/

#pragma once

class ChromChars :
	public Chromosome
{
public:
	/** chromosom - wektor gen�w typu char  */
	std::vector<char> genes;

	/** Konstruktor tworz�cy nowy wektor gen�w o rozmiarze iGenesDim */
	ChromChars(int iGenesDim);

	/** Zwraca wektor zawieraj�cy geny 	 */
	std::vector<char> getGenes();

	/** Zwraca geny w postaci string	 */
	std::string getGenesAsStr();

	/** Nadpisuje geny obecnego chromosomu poprzez kopiowania ich ze wzkazanego */
	void copyChromGenes(Chromosome* chromosome);

	/** Zwraca pojedynczy gen ze wzkazanej pozycji */
	char getGene(int iGene);

	/** Ustawia pul� gen�w na podstawie przekazanych danych w postaci string */
	void setGenesFromStr(const std::string &sChromosome);

	/** Oblicza ile gen�w jest wsp�lnych dla danego chromosomu wzgl�dem przekazanego w argumencie */
	int getNumGenesInCommon(Chromosome* chromosome);

	~ChromChars();
};

