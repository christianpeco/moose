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

#include "ElasticMaterial2DFrac.h"

template<>
InputParameters validParams<ElasticMaterial2DFrac>()
{
  InputParameters params = validParams<Material>();

 //Description
  params.addClassDescription("Elastic Material for fracture working in pure 2D. A tensile-compressive split is performed and a damage variable is included to affect the elastic constants.");

  //elastic properties
  params.addParam<Real>("bulk_modulus", "The bulk modulus for the material.");
  params.addParam<Real>("lambda", "Lame's first constant for the material.");
  params.addParam<Real>("poissons_ratio", "Poisson's ratio for the material.");
  params.addParam<Real>("shear_modulus", "The shear modulus of the material.");
  params.addParam<Real>("youngs_modulus", "Young's modulus of the material.");
  params.addParam<MaterialPropertyName>("compfactor_name", "fcomp", "A factor affecting the compressive strength");
  params.addParam<MaterialPropertyName>("tensfactor_name", "ftens", "A factor affecting the tensile strength");

  //options
  params.addParam<bool>("irrev", true, "A strain history is kept and used to compute the energies for damage");
  params.addParam<bool>("split", true, "A tensile/compressive split is performed");
  params.addParam<bool>("staggered",true, "Staggered approach. Stress and Jacobian use a staggered elastic constants");
  params.addParam<Real>("eta", "numerical parameter to prevent zero in stiffness matrix");
  params.addRequiredParam<std::string>("base_name", "Material name");

  //coupled variables
  params.addCoupledVar("damage", "Coupled variable affecting the moduli");


  return params;
}

ElasticMaterial2DFrac::ElasticMaterial2DFrac(const InputParameters & parameters) :
   DerivativeMaterialInterface<Material>(parameters),

   //parameters
   _bulk_modulus_set( parameters.isParamValid("bulk_modulus") ),
   _lambda_set( parameters.isParamValid("lambda") ),
   _poissons_ratio_set( parameters.isParamValid("poissons_ratio") ),
   _shear_modulus_set( parameters.isParamValid("shear_modulus") ),
   _youngs_modulus_set( parameters.isParamValid("youngs_modulus") ),

   _irrev_set( parameters.isParamValid("irrev") ),
   _split_set( parameters.isParamValid("split") ),
   _staggered_set( parameters.isParamValid("staggered") ),
   _eta_set( parameters.isParamValid("eta") ),

   _bulk_modulus( _bulk_modulus_set ? getParam<Real>("bulk_modulus") : -1 ),
   _lambda( _lambda_set ? getParam<Real>("lambda") : -1 ),
   _poissons_ratio( _poissons_ratio_set ?  getParam<Real>("poissons_ratio") : -1 ),
   _shear_modulus( _shear_modulus_set ? getParam<Real>("shear_modulus") : -1 ),
   _youngs_modulus( _youngs_modulus_set ? getParam<Real>("youngs_modulus") : -1 ),
   _compfactor(getMaterialProperty<Real>("compfactor_name")),
   _tensfactor(getMaterialProperty<Real>("tensfactor_name")),


   _irrev( _irrev_set ? getParam<bool>("irrev") : true),
   _split( _split_set ? getParam<bool>("split") : true),
   _staggered( _staggered_set ? getParam<bool>("staggered") : true),
   _eta( _eta_set ? getParam<Real>("eta") : -1 ),
   _base_name(getParam<std::string>("base_name")),


   _damage(coupledValue("damage")),



   //fetch strains from regular material strain
   _total_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "_" + "total_strain")),

   //declare stress and derivatives for elastic kernels
   _stress(declareProperty<std::vector<Real> >(_base_name + "_stress")),
   _Jacobian_mult(declareProperty<std::vector<Real> >(_base_name + "_Jacobian_mult")),

   //declare stateful properties for staggered approach (trace, strain plus, strain minus)
   _trace_st(declareProperty<Real>(_base_name + "_" + "trace_st")),
   _trace_st_old(declarePropertyOld<Real>(_base_name + "_" + "trace_st")),
   _strain_plus_st(declareProperty<std::vector<Real> >(_base_name + "_" + "strain_plus_st")),
   _strain_plus_st_old(declarePropertyOld<std::vector<Real> >(_base_name + "_" + "strain_plus_st")),
   _strain_minus_st(declareProperty<std::vector<Real> >(_base_name + "_" + "strain_minus_st")),
   _strain_minus_st_old(declarePropertyOld<std::vector<Real> >(_base_name + "_" + "strain_minus_st")),
   _damage_st(declareProperty<Real>(_base_name + "_" + "damage_st")),
   _damage_st_old(declarePropertyOld<Real>(_base_name + "_" + "damage_st")),

   //declare energies tensile
   _energy_tensile_raw(declareProperty<Real>(_base_name + "_" +"energy_tensile_raw")),
   _energy_tensile_raw_history(declareProperty<Real>(_base_name + "_" + "energy_tensile_raw_history")),
   _energy_tensile_raw_history_old(declarePropertyOld<Real>(_base_name + "_" + "energy_tensile_raw_history")),

   //declare energy derivatives w.r.t. damage
   _energy_tensile(declareProperty<Real>(_base_name + "_" + "energy_tensile")),
   _d_energy_tensile_dc(declareProperty<Real>(_base_name + "_" + "d_energy_tensile_dc")),
   _d2_energy_tensile_dc2(declareProperty<Real>(_base_name + "_" + "d2_energy_tensile_dc2")),

   //declare energies compressive
   _energy_compressive_raw(declareProperty<Real>(_base_name + "_" + "energy_compressive_raw")),
   _energy_compressive(declareProperty<Real>(_base_name + "_" + "energy_compressive")),

_trace_trial(declareProperty<Real>(_base_name + "_" + "trace_trial"))
{

    iso_const.resize(2);

    if (_lambda_set && _shear_modulus_set)
    {
      iso_const[0] = _lambda;
      iso_const[1] = _shear_modulus;
    }
    else if (_youngs_modulus_set && _poissons_ratio_set)
    {
      iso_const[0] = _youngs_modulus*_poissons_ratio/((1 + _poissons_ratio)*(1 - _poissons_ratio));//expression for pure 2D.
      iso_const[1] = _youngs_modulus/(2*(1 + _poissons_ratio));
    }
    else if (_shear_modulus_set && _bulk_modulus_set)
    {
      iso_const[0] = _bulk_modulus - 2.0/3.0*_shear_modulus;
      iso_const[1] = _shear_modulus;
    }
    else
      mooseError("Incorrect combination of elastic properties. Possible combinations are: lambda and shear_modulus, youngs_modulus and poissons_ratio, or bulk_modulus and shear_modulus.");

}


void
ElasticMaterial2DFrac::initQpStatefulProperties()
{
    //dimensionate the tensors (strain plus, strain minus)
    _strain_plus_st[_qp].resize(4);
    _strain_minus_st[_qp].resize(4);
    _strain_plus_st_old[_qp].resize(4);
    _strain_minus_st_old[_qp].resize(4);

    //initialize
    _trace_st[_qp]=0.0;  _damage_st[_qp]=0.0;  _energy_tensile_raw_history[_qp]=0.0;
    for (int i=0;i<4;i++){_strain_plus_st[_qp][i]=0.0;_strain_minus_st[_qp][i]=0.0;}
}


void
ElasticMaterial2DFrac::computeQpProperties()
{

    int i, j, k, l, alpha, beta;

    //identity tensor
    double iden[4]; iden[0]=1.0;iden[1]=0.0;iden[2]=0.0;iden[3]=1.0;

    //Fetch strain 2D for this qp
    int dim = 2;
    std::vector<Real> strain2d; strain2d.resize(dim*dim);
    std::vector<Real> stress2d; stress2d.resize(dim*dim);

    for (i=0;i<dim;i++){
        for (j=0;j<dim;j++){
            strain2d[i*dim+j]=_total_strain[_qp](i,j);
        }
    }

    //variables for split

    //dim values,eigenvectors and derivative of eigenvectors w.r.t. tensor entries
    Real eigvals[2];
    Real eigvecs[2][2];
    Real deigvecs[2][2][2];
    //dim and ini the the tensile(plus) and compressive(minus) parts for this qp
    Real strain_plus[dim*dim];
    Real strain_minus[dim*dim];
    // create elastic properties for tensile and compressive
    double iso_const_plus[2];
    double iso_const_minus[2];
    double trace, trace_plus, trace_minus;
     double _conc_factor, d_conc_factor_dc, d2_conc_factor_dc2, aux;
    double dstrain_plus[dim][dim][dim][dim];
    double dstrain_minus[dim][dim][dim][dim];
    //create identity4 tensor and open identity * identity, and compute _Jacob2d
    double identity4[dim*dim*dim*dim];for (i=0;i<dim*dim*dim*dim;i++) identity4[i]=0.0;
    double ideniden[dim*dim*dim*dim];for (i=0;i<dim*dim*dim*dim;i++) ideniden[i]=0.0;
    double _Jacob2d_plus[dim*dim*dim*dim];for (i=0;i<dim*dim*dim*dim;i++) _Jacob2d_plus[i]=0.0;
    double _Jacob2d_minus[dim*dim*dim*dim];for (i=0;i<dim*dim*dim*dim;i++) _Jacob2d_minus[i]=0.0;
    int pos;

    double _Jacob2d[dim*dim*dim*dim];for (i=0;i<dim*dim*dim*dim;i++) _Jacob2d[i]=0.0;



    if (_split==false)
    {

        for (i=0;i<dim*dim;i++){
            stress2d[i] = iso_const[0]*(trace*iden[i]) + 2.0*iso_const[1]*strain2d[i];
            _stress[_qp][i]=stress2d[i];
        }

        for (i=0;i<dim;i++){
            for (j=0;j<dim;j++){
                for (k=0;k<dim;k++){
                    for (l=0;l<dim;l++){

                        pos = i*(dim*dim*dim)+j*(dim*dim)+k*(dim)+l;

                        if ((i==k)&&(j==l)) identity4[pos] += 0.5;
                        if ((i==l)&&(j==k)) identity4[pos] += 0.5;

                        ideniden[pos] = iden[i*dim+j]*iden[k*dim+l];
                        _Jacob2d[pos] = iso_const[0]*(ideniden[pos]) + 2.0*iso_const[1]*identity4[pos];
                        _Jacobian_mult[_qp][pos] = _Jacob2d[pos];
                    }
                }
            }
        }

    }

 if (_split==true)
 {

 //************* STARTING TENSILE-COMPRESSIVE SPLIT (BASED ON EIGENVALUE DECOMP) **************//



    // 1. EIGVAL AND EIGVECS COMPUTATION
    // a is an auxiliar destroyed by dsyev, and if calculation_type == "V" then eigenvectors are placed there
    // Note the explicit symmeterisation
    std::vector<PetscScalar> a;
    a.resize(dim*dim);
    int lwork = 66 * dim;
    int info;
    std::vector<PetscScalar> work(lwork);
    for (i = 0; i < dim; ++i) for (j = 0; j < dim; ++j) a[i*dim + j] = 0.5 * (strain2d[i*dim + j]+ strain2d[j*dim + i]);
    // compute the eigenvalues only (if calculation_type == "N"),
    // or both the eigenvalues and eigenvectors (if calculation_type == "V")
    // assume upper triangle of a is stored (second "U")
    LAPACKsyev_("V", "U", &dim, &a[0], &dim, &eigvals[0], &work[0], &lwork, &info);
    for (i = 0; i < dim; ++i) for (j = 0; j < dim; ++j) eigvecs[j][i] = a[i*dim + j];//store eigvecs as columns.

    // 2. EIGVAL first derivatives w.r.t. tensor values (=open products of eigvecs).

    for (i=0;i<dim;i++){
        for (j=0;j<dim;j++){
            for (k=0;k<dim;k++){
                deigvecs[i][j][k] = eigvecs[j][i]*eigvecs[k][i];
            }
        }
    }

    //in case equal eigenvalues, for discontinuity
    if (eigvals[0] == eigvals[1]){
        for (i=0;i<dim;i++){
            for (j=0;j<dim;j++){
                deigvecs[0][i][j] = deigvecs[1][i][j] = (deigvecs[0][i][j] + deigvecs[1][i][j]) / 2.0;
            }
        }
    }

     // 3. USE EIGVALS AND EIGVECS TO PERFORM BASIC DECOMPOSITION

    for (i=0;i<dim*dim;i++) {strain_plus[i]=0.0;strain_minus[i]=0.0;}

    for (i=0;i<dim;i++){
        for (j=0;j<dim;j++){
            for (k=0;k<dim;k++){
                if (eigvals[i]>=0.0) strain_plus[j*dim+k] += eigvals[i]*deigvecs[i][j][k];
                if (eigvals[i]<0.0)  strain_minus[j*dim+k] += eigvals[i]*deigvecs[i][j][k];
            }
        }
    }

     // 4. COMPUTE THE STRESS PLUS AND MINUS

    for (i=0;i<dim;i++) {iso_const_plus[i] = iso_const_minus[i] = iso_const[i];}

    //tensile and compressive character of trace
    trace_plus=0.0;trace_minus=0.0;
    trace = strain2d[0] + strain2d[3];
    //if trace is positive, eliminate the corresponding part in the compressive counterpart (and opposite)
    if (trace>=0.0) {trace_plus = trace; iso_const_minus[0]=0.0;}
    if (trace<0.0) {trace_minus = trace; iso_const_plus[0]=0.0;}
    //dim stress plus and minus
    std::vector<Real> stress2d_plus; stress2d_plus.resize(dim*dim);
    std::vector<Real> stress2d_minus; stress2d_minus.resize(dim*dim);
    _stress[_qp].resize(dim*dim);

    //store elastic values to compute  elastic energy for surfactant balance with stateful strains
    _trace_st[_qp]=trace;
    for (i=0;i<dim*dim;i++) {_strain_plus_st[_qp][i]=strain_plus[i];_strain_minus_st[_qp][i]=strain_minus[i];}

    //damage tensile properties

    _damage_st[_qp] = _damage[_qp];

    if (_staggered==true) //tensile elastic variables affected by damage (last step solution)
    {
    _conc_factor = (1.0-_damage_st_old[_qp])*(1.0-_damage_st_old[_qp])+_eta;
    }
    if (_staggered==false) //tensile elastic variables affected by damage (current step solution)
    {
    _conc_factor = (1.0-_damage_st[_qp])*(1.0-_damage_st[_qp])+_eta;
    }

    for (i=0;i<dim;i++) iso_const_plus[i] = _conc_factor*_tensfactor[_qp]*iso_const_plus[i];
    for (i=0;i<dim;i++) iso_const_minus[i] = _compfactor[_qp]*iso_const_minus[i];

    //compute the stress
    for (i=0;i<dim*dim;i++){
        stress2d_plus[i] = iso_const_plus[0]*(trace_plus*iden[i]) + 2.0*iso_const_plus[1]*strain_plus[i];
        stress2d_minus[i] = iso_const_minus[0]*(trace_minus*iden[i]) + 2.0*iso_const_minus[1]*strain_minus[i];
        _stress[_qp][i]=stress2d_plus[i]+stress2d_minus[i];
    }

    // 5. COMPUTE THE JACOBIAN PLUS AND MINUS

    //Compute the second derivatives of the eigenvalues w.r.t. strain tensor
    double d2eigvals[dim][dim][dim][dim][dim];
    for (alpha = 0; alpha < dim; ++alpha){
        for (i = 0; i < dim; ++i){
            for (j = 0; j < dim; ++j){
                for (k = 0; k < dim; ++k){
                    for (l = 0; l < dim; ++l){ d2eigvals[alpha][i][j][k][l]=0.0;}}}}}

    for (alpha = 0; alpha < dim; ++alpha)
        for (beta = 0; beta < dim; ++beta)
        {
            if (eigvals[alpha] == eigvals[beta])
                continue;

            for (i = 0; i < dim; ++i)
                for (j = 0; j < dim; ++j)
                    for (k = 0; k < dim; ++k)
                        for (l = 0; l < dim; ++l)
                        {
                            d2eigvals[alpha][i][j][k][l] += 0.5 * (eigvecs[beta][i] * eigvecs[alpha][j] + eigvecs[alpha][i] * eigvecs[beta][j])
                                    * (eigvecs[beta][k] * eigvecs[alpha][l] + eigvecs[beta][l] * eigvecs[alpha][k])
                                    / (eigvals[alpha] - eigvals[beta]);
                        }
        }

    //Compute the derivative of strain_plus and strain_minus w.r.t. strain tensor

    for (i = 0; i < dim; ++i){
        for (j = 0; j < dim; ++j){
            for (k = 0; k < dim; ++k){
                for (l = 0; l < dim; ++l){ dstrain_plus[i][j][k][l]=0.0;dstrain_minus[i][j][k][l]=0.0;}}}}

    for (alpha = 0; alpha < dim; ++alpha)//loop in eigenvalues
    {
        for (i = 0; i < dim; ++i){
            for (j = 0; j < dim; ++j){
                for (k = 0; k < dim; ++k){
                    for (l = 0; l < dim; ++l){
                        if (eigvals[alpha]>=0.0) dstrain_plus[i][j][k][l] += eigvals[alpha]*d2eigvals[alpha][i][j][k][l] + deigvecs[alpha][i][j]*deigvecs[alpha][k][l];
                        if (eigvals[alpha]<0.0) dstrain_minus[i][j][k][l] += eigvals[alpha]*d2eigvals[alpha][i][j][k][l] + deigvecs[alpha][i][j]*deigvecs[alpha][k][l];
                    }
                }
            }
        }
    }


//Compute stress derivative w.r.t. strain tensor

_Jacobian_mult[_qp].resize(dim*dim*dim*dim);

for (i=0;i<dim;i++){
    for (j=0;j<dim;j++){
        for (k=0;k<dim;k++){
            for (l=0;l<dim;l++){

                pos = i*(dim*dim*dim)+j*(dim*dim)+k*(dim)+l;

                ideniden[pos] = iden[i*dim+j]*iden[k*dim+l];//CHECK!


                _Jacob2d_plus[pos] = iso_const_plus[0]*(ideniden[pos]) + 2.0*iso_const_plus[1]*dstrain_plus[i][j][k][l];
                _Jacob2d_minus[pos] = iso_const_minus[0]*(ideniden[pos]) + 2.0*iso_const_minus[1]*dstrain_minus[i][j][k][l];
                _Jacobian_mult[_qp][pos] = _Jacob2d_plus[pos]+_Jacob2d_minus[pos];
            }
        }
    }
}

 }

 //************* FINISHED TENSILE-COMPRESSIVE SPLIT (BASED ON EIGENVALUE DECOMP) **************//



 //COMPUTE TENSILE ENERGIES AND DERIVATIVES W.R.T. DAMAGE (for typical usage with phasefield phase kernel)


 //recompute tensile and compressive const to use in energy
 for (i=0;i<dim;i++) {iso_const_plus[i] = iso_const[i];iso_const_minus[i] = iso_const[i];}

 if (_staggered==true)
 {
     if (_trace_st_old[_qp]>=0.0) {trace_plus = _trace_st_old[_qp];iso_const_minus[0]=0.0;}
     if (_trace_st_old[_qp]<0.0) {trace_minus = _trace_st_old[_qp];iso_const_plus[0]=0.0;}
 }

 if (_staggered==false)
 {
     if (_trace_st[_qp]>=0.0) {trace_plus = _trace_st[_qp];iso_const_minus[0]=0.0;}
     if (_trace_st[_qp]<0.0) {trace_minus = _trace_st[_qp];iso_const_plus[0]=0.0;}
 }

 //recompute concentration factor and derivatives for surfactant damage evolution (always at current step)

 _conc_factor = (1.0-_damage[_qp])*(1.0-_damage[_qp])+_eta;
 d_conc_factor_dc = -2.0*(1.0-_damage[_qp]);
 d2_conc_factor_dc2 = 2.0;

 //TENSILE ( relevant part for fracture)

 if (_staggered==true)
 {
     aux=0.0;
     for (i=0;i<dim*dim;i++) aux += _strain_plus_st_old[_qp][i]*_strain_plus_st_old[_qp][i];
 }
 if (_staggered==false)
 {
     aux=0.0;
     for (i=0;i<dim*dim;i++) aux += _strain_plus_st[_qp][i]*_strain_plus_st[_qp][i];
 }


 //energy tensile raw (without the damaging factor)
 _energy_tensile_raw[_qp] = (0.5*iso_const_plus[0]*(trace_plus*trace_plus) + 2.0*iso_const_plus[1]*aux);


//irreversibility (H field keeping the maximum)
 if (_irrev==true)
 {
     if (_energy_tensile_raw[_qp] >= _energy_tensile_raw_history_old[_qp])
     {_energy_tensile_raw_history[_qp]=_energy_tensile_raw[_qp];}
     else {_energy_tensile_raw_history[_qp]=_energy_tensile_raw_history_old[_qp];}
 }
  if (_irrev==false) {_energy_tensile_raw_history[_qp]=_energy_tensile_raw[_qp];}


 //compute energies and derivatives for phasefield fracture
 _energy_tensile[_qp] = _conc_factor*_tensfactor[_qp]*_energy_tensile_raw_history[_qp];
 _d_energy_tensile_dc[_qp] = d_conc_factor_dc*_tensfactor[_qp]*_energy_tensile_raw_history[_qp];
 _d2_energy_tensile_dc2[_qp] = d2_conc_factor_dc2*_tensfactor[_qp]*_energy_tensile_raw_history[_qp];



 //energy compressive raw (no factor)

 if (_staggered==true)
 {
     aux=0.0;
     for (i=0;i<dim*dim;i++) aux += _strain_minus_st_old[_qp][i]*_strain_minus_st_old[_qp][i];
 }

 if (_staggered==false)
 {
     aux=0.0;
     for (i=0;i<dim*dim;i++) aux += strain_minus[i]*strain_minus[i];
 }

 _energy_compressive_raw[_qp] = (0.5*iso_const_minus[0]*(trace_minus*trace_minus) + 2.0*iso_const_minus[1]*aux);
 _energy_compressive[_qp] = _compfactor[_qp]*_energy_compressive_raw[_qp];


}



