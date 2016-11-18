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

#include "PhiUpdateAux.h"

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<PhiUpdateAux>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("update of field (phi) using the last step trace from an IC");
  params.addRequiredParam<std::string>("mat_name", "Praft material name");

  //coupled variables
  params.addCoupledVar("initial", "Initial field used in evolution");

  return params;
}

PhiUpdateAux::PhiUpdateAux(const InputParameters & parameters) :
  // You must call the constructor of the base class first
  Kernel(parameters),
  _mat_name(getParam<std::string>("mat_name")),
  _initial_phi(coupledValue("initial")),
  _trace_old(getMaterialPropertyOldByName<Real>(_mat_name + "_" + "trace_st"))

{}

Real PhiUpdateAux::computeQpResidual()
{

  return _test[_i][_qp]*(_u[_qp]-_initial_phi[_qp]/(1.0 + _trace_old[_qp]));
}

Real PhiUpdateAux::computeQpJacobian()
{
  return _test[_i][_qp]*_phi[_j][_qp];
}
