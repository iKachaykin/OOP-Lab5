#pragma once
#define EPS 0.00000001
#include <iostream>

class InterfaceVectorSpaceElement {
public:
	virtual double cube_norm() const = 0 ;
	virtual double oct_norm() const = 0;
	virtual double Euclid_norm() const = 0;
	virtual std::string to_string() const = 0;
};

std::ostream& operator<<(std::ostream&, const InterfaceVectorSpaceElement&);
