/**
Klasa Crossover definiuje dostêpne typy krzy¿owania
 */
#pragma once
class Crossover
{
public:
	/** krzy¿owanie jednopunktowe */
	static constexpr int ctOnePoint = 0;
	/** krzy¿owania dwupunktowe */
	static constexpr int ctTwoPoint = 1;
	/** krzy¿owanie jednolite */
	static constexpr int ctUniform = 2;
	/** krzy¿owanie losowe (jedno z powy¿szych) */
	static constexpr int ctRoulette = 3;
	Crossover() {};
	~Crossover() {};
};

