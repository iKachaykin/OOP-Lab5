#include "InterfaceLinearSquareSystem.h"

std::ostream & operator<<(std::ostream &os, InterfaceLinearSquareSystem &sys) {
	return os << sys.to_string().c_str();
}
