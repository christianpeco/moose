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

#ifndef ElasticMaterial2DFrac_H
#define ElasticMaterial2DFrac_H

#include "Material.h"
#include "RankTwoTensor.h"



//Forward Declarations
class ElasticMaterial2DFrac;

template<>
InputParameters validParams<ElasticMaterial2DFrac>();

/**
 * Praft material class that defines a few properties.
 */
class ElasticMaterial2DFrac : public DerivativeMaterialInterface<Material>
{
public:
  ElasticMaterial2DFrac(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void initQpStatefulProperties();


private:


  /// Elastic constants
  bool _bulk_modulus_set;
  bool _lambda_set;
  bool _poissons_ratio_set;
  bool _shear_modulus_set;
  bool _youngs_modulus_set;

  bool _irrev_set;
  bool _split_set;
  bool _staggered_set;
  bool _eta_set;

  Real _bulk_modulus;
  Real _lambda;
  Real _poissons_ratio;
  Real _shear_modulus;
  Real _youngs_modulus;

  const MaterialProperty<Real> & _compfactor;
  const MaterialProperty<Real> & _tensfactor;

  bool _irrev;
  bool _split;
  bool _staggered;
  Real _eta;

  std::string _base_name;

  //Coupled variable
  const VariableValue & _damage;

  const MaterialProperty<RankTwoTensor> & _total_strain;
  MaterialProperty<std::vector<Real> > & _stress;
  MaterialProperty<std::vector<Real> > & _Jacobian_mult;

  MaterialProperty<Real> & _trace_st;
  MaterialProperty<Real> & _trace_st_old;
  MaterialProperty<std::vector<Real> > & _strain_plus_st;
  MaterialProperty<std::vector<Real> > & _strain_plus_st_old;
  MaterialProperty<std::vector<Real> > & _strain_minus_st;
  MaterialProperty<std::vector<Real> > & _strain_minus_st_old;
  MaterialProperty<Real> & _damage_st;
  MaterialProperty<Real> & _damage_st_old;

  MaterialProperty<Real> & _energy_tensile_raw;
  MaterialProperty<Real> & _energy_tensile_raw_history;
  MaterialProperty<Real> & _energy_tensile_raw_history_old;

  MaterialProperty<Real> & _energy_tensile;
  MaterialProperty<Real> & _d_energy_tensile_dc;
  MaterialProperty<Real> & _d2_energy_tensile_dc2;

  MaterialProperty<Real> & _energy_compressive_raw;
  MaterialProperty<Real> & _energy_compressive;


  std::vector<Real> iso_const;


  MaterialProperty<Real> & _trace_trial;

};

#endif //ElasticMaterial2DFrac_H
