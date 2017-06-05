#pragma once
#include "MathVector.h"
class InterfaceLinearSquareSystem {
public:
	virtual std::string to_string() = 0;
	virtual int dimension() const = 0;
	virtual MathVector solution() const = 0;
	virtual double error() const = 0;
};

std::ostream& operator<<(std::ostream &, InterfaceLinearSquareSystem &);