/**
* Klasa rozszerza Chromosome dla chromosomów typu Char
*/

#pragma once

class ChromChars :
	public Chromosome
{
public:
	/** chromosom - wektor genów typu char  */
	std::vector<char> genes;

	/** Konstruktor tworz¹cy nowy wektor genów o rozmiarze iGenesDim */
	ChromChars(int iGenesDim);

	/** Zwraca wektor zawieraj¹cy geny 	 */
	std::vector<char> getGenes();

	/** Zwraca geny w postaci string	 */
	std::string getGenesAsStr();

	/** Nadpisuje geny obecnego chromosomu poprzez kopiowania ich ze wzkazanego */
	void copyChromGenes(Chromosome* chromosome);

	/** Zwraca pojedynczy gen ze wzkazanej pozycji */
	char getGene(int iGene);

	/** Ustawia pulê genów na podstawie przekazanych danych w postaci string */
	void setGenesFromStr(const std::string &sChromosome);

	/** Oblicza ile genów jest wspólnych dla danego chromosomu wzglêdem przekazanego w argumencie */
	int getNumGenesInCommon(Chromosome* chromosome);

	~ChromChars();
};

