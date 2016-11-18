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

#ifndef SurfactantPressure2DFrac_H
#define SurfactantPressure2DFrac_H

#include "Material.h"
#include "RankTwoTensor.h"


//Forward Declarations
class SurfactantPressure2DFrac;
class Function;

template<>
InputParameters validParams<SurfactantPressure2DFrac>();

/**
 * Praft material class that defines a few properties.
 */
class SurfactantPressure2DFrac : public DerivativeMaterialInterface<Material>
{
public:
  SurfactantPressure2DFrac(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void initQpStatefulProperties();

private:
  /**
   * parameters
  */

  bool _const_set;
  bool  _pres_type_set;
  bool _irrev_set;
  bool _staggered_set;

  Real _const;
  int _pres_type;
  bool _irrev;
  bool _staggered;

  Function & _loadingfunction;
  std::string _base_name;


  /**
   * This is the member reference that will hold the computed values
   * for the values property and gradients in this class.
   */



  /**
   * This is the member reference that will hold the gradient
   * of the coupled variable
   */



  const VariableValue & _concentration;
  const VariableValue & _concentration_ini;
  const MaterialProperty<Real> & _trace;
  const MaterialProperty<Real> & _trace_old;

  MaterialProperty<Real> &_concentration_history;
  MaterialProperty<Real> &_concentration_history_old;
  MaterialProperty<Real> &_trace_history;
  MaterialProperty<Real> &_trace_history_old;


  MaterialProperty<Real> & _overisostress;
  MaterialProperty<Real> & _energy_over;
  MaterialProperty<Real> & _d_energy_over_dc;
  MaterialProperty<Real> & _d2_energy_over_dc2;





};

#endif //SurfactantPressure2DFrac_H
