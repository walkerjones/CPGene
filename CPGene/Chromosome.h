/**
 * Chromosome to abstrakcyjna klasa bazowa dla wszystkich typów chromosomów. It defines each chromosome's
 * Dla ka¿dego chromosomu deklaruje ona geny, dopasowanie, pozycjê w rankingu dopasowania,
 * oraz proste metody kopiowania i zwracania wartoœci chromosomów w psotaci string
 */
#pragma once
#include <string>
class Chromosome
{
public:
	double fitness = 0.0;
	int fitnessRank = 0;
	virtual std::string getGenesAsStr()=0;
	virtual void copyChromGenes(Chromosome *chromosome)=0;
	virtual int getNumGenesInCommon(Chromosome *chromosome)=0;
	Chromosome() {};
	~Chromosome() {};
};

