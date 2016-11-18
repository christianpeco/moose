/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ACInterfaceFrac_H
#define ACInterfaceFrac_H

#include "Kernel.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

class ACInterfaceFrac;

template<>
InputParameters validParams<ACInterfaceFrac>();

/**
 * Compute the Allen-Cahn interface term with the weak form residual
 * \f$ \left( \kappa_i \nabla\eta_i, \nabla (L_i \psi) \right) \f$
 */
class ACInterfaceFrac : public DerivativeMaterialInterface<JvarMapKernelInterface<Kernel> >
{
public:
  ACInterfaceFrac(const InputParameters & parameters);
  virtual void initialSetup();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  RealGradient gradL();
  RealGradient gradKappa();
  RealGradient gradGc();

  /// the \f$ \kappa\nabla(L\psi) \f$ term
  RealGradient kappaNablaLGcPsi();

  /// Mobility
  const MaterialProperty<Real> & _L;
  /// Interfacial parameter
  const MaterialProperty<Real> & _kappa;
  /// Toughness
  const MaterialProperty<Real> & _Gc;

  /// flag set if L is a function of non-linear variables in args
  const bool _variable_L;
  /// flag set if Gc is a function of non-linear variables in args
  const bool _variable_Gc;



  /// @{ Mobility derivatives w.r.t. order parameter
  const MaterialProperty<Real> & _dLdop;
  const MaterialProperty<Real> & _d2Ldop2;
  /// @}

  /// kappa derivative w.r.t. order parameter
  const MaterialProperty<Real> & _dkappadop;

  /// @{ Toughness derivatives w.r.t. order parameter
  const MaterialProperty<Real> & _dGcdop;
  const MaterialProperty<Real> & _d2Gcdop2;
  /// @}

  /// number of coupled variables
  const unsigned int _nvar;

  /// @{ Mobility derivative w.r.t. other coupled variables
  std::vector<const MaterialProperty<Real> *> _dLdarg;
  std::vector<const MaterialProperty<Real> *> _d2Ldargdop;
  std::vector<std::vector<const MaterialProperty<Real> *> > _d2Ldarg2;
  /// @}

  /// kappa derivative w.r.t. other coupled variables
  std::vector<const MaterialProperty<Real> *> _dkappadarg;

  /// @{ Toughness derivative w.r.t. other coupled variables
  std::vector<const MaterialProperty<Real> *> _dGcdarg;
  std::vector<const MaterialProperty<Real> *> _d2Gcdargdop;
  std::vector<std::vector<const MaterialProperty<Real> *> > _d2Gcdarg2;
  /// @}

  /// Gradients for all coupled variables
  std::vector<const VariableGradient *> _gradarg;
};

#endif //ACInterfaceFrac_H
