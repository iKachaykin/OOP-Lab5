#pragma once
#include "InterfaceVectorSpaceElement.h"
class InterfaceMatrix 
	: public InterfaceVectorSpaceElement {
public:
	virtual int rows() const = 0;
	virtual int cols() const = 0;
	virtual void transpose() = 0;
};

