/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is, 
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#ifndef _TRANSFORMATIONS_H_
#define _TRANSFORMATIONS_H_
#include <iostream>
#include <string>
#include "types.h"

using namespace std;

/**
  * Abstract class that represents a transformation from n bits to n bits.
  */
class Transformation {
public:
    Transformation() {};
    /**
      * Abstract method that returns the number of bits of its domain and range.
      */
    virtual unsigned int getWidth() const = 0;
    /**
      * Abstract method that applies the transformation onto the parameter 
      * @a state.
      *
      * @param  state   A buffer on which to apply the transformation.
      *                 The state must have a size of at least 
      *                 ceil(getWidth()/8.0) bytes.
      */
    virtual void operator()(UINT8 * state) const = 0;
    /**
      * Abstract method that returns a string with a description of itself.
      */
    virtual string getDescription() const = 0;
    /**
      * Method that prints a brief description of the transformation.
      */
    friend ostream& operator<<(ostream& a, const Transformation& transformation);
};

/**
  * Abstract class that represents a permutation from n bits to n bits.
  */
class Permutation : public Transformation {
public:
    Permutation() : Transformation() {};
    /**
      * Abstract method that applies the <em>inverse</em> of the permutation 
      * onto the parameter @a state.
      *
      * @param  state   A buffer on which to apply the inverse permutation.
      *                 The state must have a size of at least 
      *                 ceil(getWidth()/8.0) bytes.
      */
    virtual void inverse(UINT8 * state) const = 0;
};

/**
  * Class that implements the simplest possible permutation: the identity.
  */
class Identity : public Permutation {
protected:
    unsigned int width;
public:
    Identity(unsigned int aWidth) : Permutation(), width(aWidth) {};
    virtual unsigned int getWidth() const { return width; }
    virtual void operator()(UINT8 * state) const {}
    virtual string getDescription() const { return "Identity"; }
    virtual void inverse(UINT8 * state) const {}
};

#endif
