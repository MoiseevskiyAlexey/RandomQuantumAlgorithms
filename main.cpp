#include <iostream>
#include <QuEST.h>
#include "randqalg.h"

#include <vector>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <cmath>

#define LINES 3
#define COLS 4
#define DEPTH 40
#define AVG 200
#define ALGS 50

using namespace std;

int main()
{
	//freopen("AllErrsExtreme2.txt", "w", stdout);
	srand( (unsigned)time(NULL) );
	QubitArray qubits;

	qubits.resize(COLS - 1, LINES -1);
	qubits.resize(COLS, LINES);

	qubits.setSingleErrRate(0.001);
	qubits.setMultiErrRate(0.01);
	qubits.setEnvCoupling(0.001);
	qubits.setSingleGateTime(0.1);
	qubits.setMultiGateTime(1.0);
	qubits.setLoseTime(10000);
	qubits.setDynamicNoise(0.04);
	qubits.setSpamErr0to1(0.5);
	qubits.setSpamErr1to0(0.5);

	RandQAlg alg(COLS, LINES, DEPTH);

	constexpr int ampNum = (1 << LINES * COLS);
	std::vector<double> amps(ampNum, 0.);

	cout << fixed << "{ ";
	for(int k = 0; k < ALGS; k++)
	{
		cerr << "Step " << k + 1 << " of " << ALGS << "\n";
		for(int i = 0; i < ampNum; i++)
			amps[i] = 0.0;
		alg.generate();
		for(int i = 0; i < AVG; i++)
		{
			alg.evaluate(qubits);
			for(int j = 0; j < ampNum; j++)
				amps[j] += qubits.getSquaredAmp(j);
		}

		for(int j = 0; j < ampNum; j++)
		{
			cout << amps[j] * ampNum / AVG;
			if(j + 1 < ampNum || k + 1 < ALGS)
				cout << ",\n";
		}
	}
	cout << "}";

	return 0;
}
