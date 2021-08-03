#ifndef QUBITARRAY_H
#define QUBITARRAY_H

#include <QuEST.h>
#include <random>
#include <vector>
#include <functional>

#define _USE_MATH_DEFINES
#define OUTCOME_PROB_EPS 1e-10

struct cords
{
	unsigned x;
	unsigned y;
};

class QubitArray
{
public:
	QubitArray(unsigned a = 1, unsigned b = 1);
	~QubitArray();

	void setEnvCoupling(double val);
	void setSingleGateCoupling(double val);
	void setMultiGateCoupling(double val);
	void setSingleErrRate(double val);
	void setMultiErrRate(double val);
	void setSingleGateTime(double val);
	void setMultiGateTime(double val);
	void setLoseTime(double val);
	void setDynamicNoise(double val);
	void setSpamErr0to1(double val);
	void setSpamErr1to0(double val);
	void setSpamErr(double val){ setSpamErr0to1(val); setSpamErr1to0(val); }
	void setAmpDampingRate(double val);
	void setCZFail(double val);
	void setRydbergLoseTime(double val);

	unsigned getXSize(){ return xSize; }
	unsigned getYSize(){ return ySize; }

	void reset();
	void resize(unsigned newX, unsigned newY);
	void generateBell(cords first, cords sec);
	void pureX(cords target){ pauliX(qubits, getIndex(target)); }
	void pureH(cords target){ hadamard(qubits, getIndex(target)); }
	void pureCZ(cords control, cords target){ controlledPhaseFlip(qubits, getIndex(control), getIndex(target)); }
	void dropQubit(unsigned index);
	void cz(cords control, cords target);
	void swap(cords first, cords sec);
	int move(cords init, cords dest);

	void applySingleGate(cords target, std::function<void(Qureg, unsigned)> gate);
	void applyRotation(cords target, Vector v, double angle);
	void hadamardGate(cords target){ applySingleGate(target, hadamard); }
	void sqrtX(cords target){ applySingleGate(target, [](Qureg reg, unsigned index){rotateX(reg, index, M_PI_2);}); }
	void sqrtY(cords target){ applySingleGate(target, [](Qureg reg, unsigned index){rotateY(reg, index, M_PI_2);}); }
	void TGate(cords target){ applySingleGate(target, tGate); }
	void XGate(cords target){ applySingleGate(target, pauliX); }

	double calcBellFidelity(cords first, cords sec);
	double calcBellFidelityDirect(cords first, cords sec);
	int meas(cords target);
	double getSquaredAmp(unsigned index);
	friend double fidelity(const QubitArray &arr1, const QubitArray &arr2){ return calcFidelity(arr1.qubits, arr2.qubits); }
	double calcProb(cords target);

	unsigned getIndex(cords c) const;
	cords getCords(unsigned index) const;

private:
	QuESTEnv env;
	Qureg qubits;

	unsigned xSize;
	unsigned ySize;
	double envCoupling;
	double singleGateCoupling;
	double multiGateCoupling;
	double singleErrRate;
	double multiErrRate;
	double singleGateTime;
	double multiGateTime;
	double loseTime;
	double dynamicNoise;
	double spamError1to0;
	double spamError0to1;
	double ampDampingRate;
	double czFailRate;
	double rydbergLoseTime;

	void updateTotalTime();
	std::random_device rd{};
	std::mt19937 gen{rd()};

	void applyNoiseGate(unsigned index, double coupling, double time);
	void applyNoise(unsigned index);
	void applyNoise(cords target){ applyNoise(getIndex(target)); }
	void applyNoise();
	void applySingleGateErr(cords target);
	void applyMultiGateErr(cords target);
	void applyMultiGateErr(cords first, cords sec){ applyMultiGateErr(first); applyMultiGateErr(sec); }
	void applyDamping(unsigned index, double time);

	double totalTime;
	bool singleGateInCurLayer;
	bool multiGateInCurLayer;
	std::vector <int> lastNoiseTime;
	std::vector <bool> usedInCurLayer;
	std::vector <bool> isLost;
	std::vector <bool> isRydberg;
	void startNewLayer();
};

#endif // QUBITARRAY_H
