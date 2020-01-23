/*	Plik z przyk³adami zastosowania biblioteki.	Wybór przyk³adu dostêpny poprzez 
	odkomentowanie (tylko) jednej z poni¿szych definicji preprocesora */

//#define TrigFun 1
//#define BinaryOnes 1 
//#define CurveFit 1
#define Knapsack 1

#include "PCH.h"

#ifdef TrigFun
double GAFloat::getFitness(int iChromIndex) { return 0.0; };
double GAString::getFitness(int iChromIndex)
{
	double value = chromStrToFloat((this->getChromosome(iChromIndex))->getGenesAsStr(), this->chromDecPts);
	value = std::sin(value) + std::cos(value);
	return (value);
}

int main()
{
	GAString trig(
		7, //liczba genów w chromosomie
		50, //rozmiar N populacji chromosomów
		0.7, //prawdopodobieñstwo krzy¿owania
		10, //szansa na losow¹ selekcjê w % (niezale¿nie od funkcji fit)
		150, //limit liczby pokoleñ
		5, //liczba wstêpnych uruchomieñ
		20, //max pokoleñ we wstêpnych uruchomieniach
		0.003, //szansa mutacji chromosomu
		6, //liczba miejsc po przecinku w chromosomie
		//jeœli chromosom ma 7 genów i 6 miejsc po przecinku to liczby maj¹ postaæ "0.123456"
		"0123456789", //pula mo¿liwych genów
		Crossover::ctTwoPoint, //typ krzy¿owania
		true //obliczanie statystyk
	);
	trig.initPopulation();
	trig.run();
	return 0;
}
#endif

#ifdef BinaryOnes
double GAFloat::getFitness(int iChromIndex) { return 0.0; };
double GAString::getFitness(int iChromIndex)
{
	std::string s = this->getChromosome(iChromIndex)->getGenesAsStr();
	return (getChromValAsDouble(s));
}

int main()
{
	GAString bin(
		20, //liczba genów chromosomu w Char
		100, //rozmiar populacji N chromosomów
		0.7, //prawdopodobieñstwo krzy¿owania
		5, //szansa na losow¹ selekcjê w % (niezale¿nie od funkcji fit)
		50, //limit liczby pokoleñ
		0, //liczba wstêpnych uruchomieñ
		10, //max pokoleñ we wstêpnych uruchomieniach
		0.01, //szansa mutacji chromosomu
		0, //liczba miejsc po przecinku w chromosomie (0 oznacza liczby ca³kowite)
		"01", //pula mo¿liwych genów
		Crossover::ctTwoPoint, //typ krzy¿owania
		true); //obliczanie statystyk
	bin.initPopulation();
	bin.run();
	return 0;
}

#endif

#ifdef CurveFit
double GAString::getFitness(int iChromIndex) { return 0.0; };
double GAFloat::getFitness(int iChromIndex)
{
	double rError = 0;
	double rPower = 0, rValue = 0;

	for (int iCurvePt = 0; iCurvePt < curveDim; iCurvePt++)
	{
		rValue = 0;

		//calculate the value of the function plugging in the genes
		//in place of coefficients (C1..C6)
		for (int iGene = 0; iGene < chromosomeDim; iGene++)
		{
			rPower = std::pow(static_cast<double>(iCurvePt), static_cast<double>(chromosomeDim) - 1 - iGene);

			rValue += this->getChromosome(iChromIndex)->getGene(iGene) * rPower;
		}

		rError = rError + std::abs(curveData[iCurvePt] - rValue);
	}

	if (std::abs(rError) > 1e-12)
	{
		return (1 / rError); //this minimizes error (find smallest error)
	}
	else
	{
		return (1 / 1e-12); //prevents divide by zero error
	}
};
void GAFloat::setCurveData(std::vector<double> &CurveData)
{
	curveDim = CurveData.size();
	curveData = std::vector<double>(curveDim);

	for (int i = 0; i < curveDim; i++)
	{
		this->curveData[i] = CurveData[i];
	}
}


int main()
{
	std::vector<double> curveData = { 1.01, 0, 1, 3.98, 9.01, 16.01, 25, 35.99 };

	GAFloat cur(
		3, // liczba genów chromosomu
		100, //rozmiar populacji N chromosomów
		0.7, //prawdopodobieñstwo krzy¿owania
		6, //szansa na losow¹ selekcjê w % (niezale¿nie od funkcji fit)
		3000, //limit liczby pokoleñ
		10, //liczba wstêpnych uruchomieñ (to build good breeding stock for final--full run)
		20, //max pokoleñ we wstêpnych uruchomieniach
		0.1, //szansa mutacji chromosomu
		Crossover::ctTwoPoint, //typ krzy¿owania
		2, //liczba miejsc po przecinku precyzji
		//jeœli chromosom ma 3 geny i 2 miejsca po przecinku, wynik ma postaæ "0.12"
		false, //wy³¹cznie liczby nieujemne
		true); //obliczanie statystyk

	cur.setCurveData(curveData);
	cur.initPopulation();
	cur.run();
	return 0;
}
#endif

#ifdef Knapsack
double GAFloat::getFitness(int iChromIndex) { return 0.0; };
double GAString::getFitness(int iChromIndex)
{
	weight= 0.0, value = 0.0;
	std::string s = this->getChromosome(iChromIndex)->getGenesAsStr();
	for (size_t i = 0; i < s.size(); i++)
	{
		if (s[i] == '1')
		{
			weight += knapWeight[i];
			value += knapValue[i];
		}
	}
	if (weight <= limit) {
		return (value);
	}
	else return 0.0;
}

void GAString::printResult(std::ofstream &txtout)
{
	weight = 0.0, value = 0.0;
	std::string s = this->chromosomes[this->bestFitnessChromIndex]->getGenesAsStr();
	for (size_t i = 0; i < s.size(); i++)
	{
		if (s[i] == '1')
		{
			txtout << std::setprecision(15) << knapWeight[i] << " " << knapValue[i] << std::endl;
			weight += knapWeight[i];
			value += knapValue[i];
		}
		else
		{
			txtout << std::endl;
		}
	}
	txtout << std::endl << "Suma wag = " << weight << " Wartosc = " << value << std::endl;
}

int main()
{
	std::ifstream txtin("input.txt");
	std::ofstream txtout("output.txt");
	std::string line;
	std::vector<double> kweight, kvalue;
	double limit = 10000000.0;
	int n = 0;
	while (std::getline(txtin, line))
	{
		std::istringstream iss(line);
		double a, b;
		if (!(iss >> a >> b)) { break; }
		kweight.push_back(a);
		kvalue.push_back(b);
		n++;
	}

	GAString bag(
		n, //liczba genów chromosomu w Char
		100, //rozmiar populacji N chromosomów
		0.7, //prawdopodobieñstwo krzy¿owania
		5, //szansa na losow¹ selekcjê w % (niezale¿nie od funkcji fit)
		50, //limit liczby pokoleñ
		0, //liczba wstêpnych uruchomieñ
		10, //max pokoleñ we wstêpnych uruchomieniach
		0.01, //szansa mutacji chromosomu
		0, //liczba miejsc po przecinku w chromosomie (0 oznacza liczby ca³kowite)
		"01", //pula mo¿liwych genów
		Crossover::ctTwoPoint, //typ krzy¿owania
		true); //obliczanie statystyk

	bag.knapWeight = kweight;
	bag.knapValue = kvalue;
	bag.limit = limit;
	bag.initPopulation();
	bag.run();
	bag.printResult(txtout);

	
	return 0;
}
#endif