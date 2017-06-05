#include "InterfaceVectorSpaceElement.h"

std::ostream & operator<<(std::ostream &os, const InterfaceVectorSpaceElement &vect) {
	os << vect.to_string().c_str();
	return os;
}

