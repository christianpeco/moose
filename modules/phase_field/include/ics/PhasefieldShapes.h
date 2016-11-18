/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef PHASEFIELDSHAPES_H
#define PHASEFIELDSHAPES_H

// MOOSE Includes
#include "InitialCondition.h"

// Forward Declarations
class PhasefieldShapes;

template<>
InputParameters validParams<PhasefieldShapes>();

/**
 * PhasefieldShapes just returns a constant value.
 */
class PhasefieldShapes : public InitialCondition
{
public:

  /**
   * Constructor: Same as the rest of the MOOSE Objects
   */
  PhasefieldShapes(const InputParameters & parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p);

private:


  Real _epsilon;
  std::vector<Real> _centers_tips;
  std::vector<std::string> _shapes;
  std::vector<std::string> _gradients;
  std::vector<Real> _gradient_values;

};

#endif //PHASEFIELDSHAPES_H
