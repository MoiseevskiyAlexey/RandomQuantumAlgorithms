#include "qubitarray.h"
#include <stdexcept>
#include <cmath>
#include <string>
#include <iostream>

QubitArray::QubitArray(unsigned a, unsigned b)
{
	if(a <= 0 || b <= 0)
		throw std::invalid_argument(std::string("Invalid Qureg size: ") +
									std::to_string(a) + " * " + std::to_string(b));
	env = createQuESTEnv();
	unsigned long seed = gen();
	seedQuEST(&seed, 1);
	xSize = a;
	ySize = b;

	qubits = createQureg(xSize * ySize, env);
	lastNoiseTime.resize(xSize * ySize);
	usedInCurLayer.resize(xSize * ySize);
	isLost.resize(xSize * ySize);
	reset();

	envCoupling = 0.2;
	singleGateCoupling = 0.2;
	multiGateCoupling = 0.2;

	setSingleErrRate(0.0);
	setMultiErrRate(0.0);
	setSpamErr(0.0);

	singleGateTime = 0.00;
	multiGateTime = 0.00;

	loseTime = 0;
	dynamicNoise = 0;
	ampDampingRate = 0;
}

QubitArray::~QubitArray()
{
	destroyQureg(qubits, env);
	destroyQuESTEnv(env);
}

unsigned QubitArray::getIndex(cords c) const
{
	if(c.x >= xSize || c.y >= ySize)
		throw std::out_of_range("Out of qubit register (x = " + std::to_string(c.x) +
								", y = " + std::to_string(c.y) + ")");
	return c.y * xSize + c.x;
}

cords QubitArray::getCords(unsigned index) const
{
	if(index >= xSize * ySize)
		throw std::out_of_range("Out of qubit register (index = " + std::to_string(index) + ")");
	return {index % xSize, index / xSize};
}

double QubitArray::getSquaredAmp(unsigned index)
{
	//supports SPAM errors
	unsigned n = xSize * ySize, spamedIndex = 0;
	double normRatio = 1.;
	std::uniform_real_distribution rand(0.0, 1.0);
	for(unsigned i = 0; i < n; i++)
	{
		bool bit = index & (1 << i);
		double r = bit ? (1. + spamError0to1 - spamError1to0) : (1. - spamError0to1 + spamError1to0);
		double errSeed = rand(gen) * r;
		normRatio *= r;
		if(bit ? (errSeed >= spamError0to1) : (errSeed < spamError1to0))
			spamedIndex += (1 << i);
	}
	double real = getAmp(qubits, spamedIndex).real, imag = getAmp(qubits, spamedIndex).imag;
	return normRatio * (real * real + imag * imag);
}

void QubitArray::dropQubit(unsigned index)
{
	if(calcProbOfOutcome(qubits, index, 0) < OUTCOME_PROB_EPS)
		pauliX(qubits, index);
	collapseToOutcome(qubits, index, 0);
}

void QubitArray::updateTotalTime()
{
	if(singleGateInCurLayer && !multiGateInCurLayer)
		totalTime += singleGateTime;
	if(multiGateInCurLayer)
		totalTime += multiGateTime;
}

void QubitArray::reset()
{
	initZeroState(qubits);
	totalTime = 0;
	std::fill(lastNoiseTime.begin(), lastNoiseTime.end(), 0);
	std::fill(isLost.begin(), isLost.end(), false);
	std::fill(usedInCurLayer.begin(), usedInCurLayer.end(), false);
	singleGateInCurLayer = false;
	multiGateInCurLayer = false;
	startNewLayer();
}

void QubitArray::resize(unsigned newX, unsigned newY)
{
	if(newX == xSize && newY == ySize)
		return;
	if(newX <= 0 || newY <= 0)
		throw std::invalid_argument(std::string("Invalid Qureg size: ") +
									std::to_string(newX) + " * " + std::to_string(newY));
	xSize = newX;
	ySize = newY;
	destroyQureg(qubits, env);
	qubits = createQureg(xSize * ySize, env);
	lastNoiseTime.resize(xSize * ySize);
	usedInCurLayer.resize(xSize * ySize);
	isLost.resize(xSize * ySize);
	reset();
}

void QubitArray::startNewLayer()
{
	updateTotalTime();
	singleGateInCurLayer = false;
	multiGateInCurLayer = false;
	std::fill(usedInCurLayer.begin(), usedInCurLayer.end(), false);
}

void QubitArray::generateBell(cords first, cords sec)
{
	hadamard(qubits, getIndex(first));
	controlledNot(qubits, getIndex(first), getIndex(sec));
	//multiGateInCurLayer = true;
	//usedInCurLayer[getIndex(first)] = usedInCurLayer[getIndex(sec)] = true;
	//lastNoiseTime[getIndex(first)] = lastNoiseTime[getIndex(sec)] = totalTime + multiGateTime;
}

void QubitArray::cz(cords control, cords target)
{
	if(getIndex(control) == getIndex(target))
		throw std::invalid_argument("In CZ instruction the control qubit cannot be equal to the target qubit");

	if(usedInCurLayer[getIndex(control)] || usedInCurLayer[getIndex(target)])
		startNewLayer();
	usedInCurLayer[getIndex(control)] = usedInCurLayer[getIndex(target)] = true;

	applyNoise(control);
	applyNoise(target);
	applyMultiGateErr(control, target);

	controlledPhaseFlip(qubits, getIndex(control), getIndex(target));

	if(isLost[getIndex(target)])
		dropQubit(getIndex(target));

	multiGateInCurLayer = true;
	applyMultiGateErr(control, target);
}

void QubitArray::swap(cords first, cords sec)
{
	//doesnt support atom loss
	if(usedInCurLayer[getIndex(first)] || usedInCurLayer[getIndex(sec)])
		startNewLayer();
	usedInCurLayer[getIndex(first)] = usedInCurLayer[getIndex(sec)] = true;

	applyNoise(first);
	applyNoise(sec);
	applyMultiGateErr(first, sec);

	swapGate(qubits, getIndex(first), getIndex(sec));

	multiGateInCurLayer = true;
	applyMultiGateErr(first, sec);
}

int QubitArray::move(cords init, cords dest)
{
	//doesnt support atom loss
	if(dest.x > xSize || dest.y > ySize || init.x > xSize || init.y > ySize)
		throw std::out_of_range("Given index is out of qubit register ranges");
	cords cur = init;
	int xStep = 1, yStep = 1;

	if (dest.x < init.x)
		xStep = -1;
	if(dest.y < init.y)
		yStep = -1;

	while (cur.y != dest.y)
	{
		swap(cur, {cur.x, cur.y + yStep});
		cur.y += yStep;
	}
	while (cur.x != dest.x)
	{
		swap(cur, {cur.x + xStep, cur.y});
		cur.x += xStep;
	}
	return abs(static_cast<int>(init.x - dest.x)) + abs(static_cast<int>(init.y - dest.y));
}

void QubitArray::applySingleGate(cords target, std::function<void(Qureg, unsigned)> gate)
{
	unsigned index = getIndex(target);
	if(usedInCurLayer[index])
		startNewLayer();
	usedInCurLayer[index] = true;
	singleGateInCurLayer = true;

	applyNoise(target);
	if(isLost[index])
		return;
	applySingleGateErr(target);
	gate(qubits, index);

	applySingleGateErr(target);
}

void QubitArray::applyRotation(cords target, Vector v, double angle)
{
	//provides dynamic noise
	applySingleGate(target,
		[angle, v](Qureg reg, unsigned index)
		{
			rotateAroundAxis(reg, index, angle, v);
		});
	std::normal_distribution<double> rand(0, dynamicNoise);
	double angleErr = rand(gen);
	rotateAroundAxis(qubits, getIndex(target), angleErr, v);
}

double QubitArray::calcBellFidelity(cords first, cords sec)
{
	updateTotalTime();
	applyNoise();
	controlledNot(qubits, getIndex(first), getIndex(sec));
	hadamard(qubits, getIndex(first));
	Complex amp = getAmp(qubits, 0);
	//hadamard(qubits, getIndex(first));
	//controlledNot(qubits, getIndex(first), getIndex(sec));
	return amp.real * amp.real + amp.imag * amp.imag;
}

double QubitArray::calcBellFidelityDirect(cords first, cords sec)
{
	updateTotalTime();
	applyNoise();

	//preparing ideal register
	Qureg noNoiseReg = createQureg(xSize * ySize, env);
	hadamard(noNoiseReg, 0);
	controlledNot(noNoiseReg, 0, 1);
	if(getIndex(first) != 0)
		swapGate(noNoiseReg, 0, getIndex(first));
	if(getIndex(sec) != 1)
		swapGate(noNoiseReg, 1, getIndex(sec));

	return calcFidelity(qubits, noNoiseReg);
}

int QubitArray::meas(cords target)
{
	int result = measure(qubits, getIndex(target));
	std::uniform_real_distribution rand(0.0, 1.0);
	double errSeed = rand(gen);
	if((result && errSeed >= spamError1to0) || ((!result) &&  errSeed < spamError0to1))
		return 1;
	if((result && errSeed < spamError1to0) || ((!result) &&  errSeed >= spamError0to1))
		return 0;
	return result;	//just to escape warning
}

double QubitArray::calcProb(cords target)
{
	double result = calcProbOfOutcome(qubits, getIndex(target), 1);
	return (result * (1 - spamError1to0)) + ((1 - result) * spamError0to1);
}

void QubitArray::applyNoiseGate(unsigned index, double coupling, double time)
{
	std::normal_distribution<double> rand(0, 1 * sqrt(time));	// \sigma = 1
	double w = rand(gen);
	rotateZ(qubits, index, 2 * sqrt(coupling) * w);
}

void QubitArray::applyNoise(unsigned index)
{
	if(isLost[index])
		return;
	std::uniform_real_distribution rand(0.0, 1.0);
	if(loseTime > 0)
		if(rand(gen) > std::exp((-1) * totalTime/loseTime))
		{
			isLost[index] = true;
			dropQubit(index);
			return;
		}
	applyNoiseGate(index, envCoupling, totalTime - lastNoiseTime[index]);
	applyDamping(index, totalTime - lastNoiseTime[index]);
	lastNoiseTime[index] = totalTime;
}

void QubitArray::applyNoise()
{
	for(unsigned i = 0; i < xSize * ySize; i++)
		applyNoise(i);
}

void QubitArray::applySingleGateErr(cords target)
{
	double effectiveTime = -(log(1 - 2 * singleErrRate) / (2 * singleGateCoupling));
	applyNoiseGate(getIndex(target), singleGateCoupling, effectiveTime);
	lastNoiseTime[getIndex(target)] = totalTime + singleGateTime;
}

void QubitArray::applyMultiGateErr(cords target)
{
	if(isLost[getIndex(target)])
		return;
	double effectiveTime = -(log(1 - 2 * multiErrRate) / (2 * multiGateCoupling));
	applyNoiseGate(getIndex(target), multiGateCoupling, effectiveTime);
	lastNoiseTime[getIndex(target)] = totalTime + multiGateTime;
}

void QubitArray::applyDamping(unsigned index, double time)
{
	std::normal_distribution<double> rand(0, 1 - exp(- ampDampingRate * time));
	double phi = rand(gen);

	ComplexMatrix2 dampingNoiseMatrix;
	dampingNoiseMatrix.real[0][0] = 1;
	dampingNoiseMatrix.real[0][1] = 0;
	dampingNoiseMatrix.real[1][0] = 0;
	dampingNoiseMatrix.real[1][1] = exp(-ampDampingRate * time / 2);

	dampingNoiseMatrix.imag[0][0] = 0;
	dampingNoiseMatrix.imag[0][1] = 0;
	dampingNoiseMatrix.imag[1][0] = phi;
	dampingNoiseMatrix.imag[1][1] = 0;

	applyMatrix2(qubits, index, dampingNoiseMatrix);
}

void QubitArray::setSingleErrRate(double val)
{
	if(val < 0 || val >= 0.5)
		throw std::invalid_argument("Single-qubit gate error rate must be in range [0, 0.5), " +
									std::to_string(val) + " obtained");
	singleErrRate = val / 2;
}

void QubitArray::setMultiErrRate(double val)
{
	if(val < 0 || val >= 0.5)
		throw std::invalid_argument("Multi-qubit gate error rate must be in range [0, 0.5), " +
									std::to_string(val) + " obtained");
	multiErrRate = val / 2;
}

void QubitArray::setSingleGateTime(double val)
{
	if(val < 0)
		throw std::invalid_argument("Single-qubit gate time must be non-negative, " +
									std::to_string(val) + " obtained");
	if (singleGateTime == val)
		return;
	singleGateTime = val;
	if(singleGateInCurLayer || multiGateInCurLayer)
		startNewLayer();
}

void QubitArray::setMultiGateTime(double val)
{

	if(val < 0)
		throw std::invalid_argument("Multi-qubit gate time must be non-negative, " +
									std::to_string(val) + " obtained");
	if(multiGateTime == val)
		return;
	multiGateTime = val;
	if(singleGateInCurLayer || multiGateInCurLayer)
		startNewLayer();
}

void QubitArray::setSingleGateCoupling(double val)
{
	if(val <= 0)
		throw std::invalid_argument("Single-qubit gate coupling constant must be positive, " +
									std::to_string(val) + " obtained");
	singleGateCoupling = val;
}

void QubitArray::setMultiGateCoupling(double val)
{
	if(val <= 0)
		throw std::invalid_argument("Multi-qubit gate coupling constant must be positive, " +
									std::to_string(val) + " obtained");
	multiGateCoupling = val;
}

void QubitArray::setEnvCoupling(double val)
{
	if(val < 0)
		throw std::invalid_argument("Environment coupling must be non-negative, " +
									std::to_string(val) + " obtained");
	envCoupling = val;
}

void QubitArray::setLoseTime(double val)
{
	if(val < 0)
		throw std::invalid_argument("Qubit life time must be non-negative, " +
									std::to_string(val) + " obtained");
	loseTime = val;
}

void QubitArray::setDynamicNoise(double val)
{
	if(val < 0)
		throw std::invalid_argument("Unitary rotation error must be non-negative, " +
									std::to_string(val) + " obtained");
	dynamicNoise = val;
}

void QubitArray::setSpamErr0to1(double val)
{
	if(val < 0 || val > 1)
		throw std::invalid_argument("SPAM error probability must be in range [0, 1], " +
									std::to_string(val) + " obtained");
	spamError0to1 = val;
}

void QubitArray::setSpamErr1to0(double val)
{
	if(val < 0 || val > 1)
		throw std::invalid_argument("SPAM error probability must be in range [0, 1], " +
									std::to_string(val) + " obtained");
	spamError1to0 = val;
}

void QubitArray::setAmpDampingRate(double val)
{
	if(val < 0 || val > 1)
		throw std::invalid_argument("Amplitude damping rate must be in range [0, 1], " +
									std::to_string(val) + " obtained");
	ampDampingRate = val;
}
