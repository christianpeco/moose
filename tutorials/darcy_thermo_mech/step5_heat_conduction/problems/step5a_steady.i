[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 10
  xmax = 0.304 # Length of test chamber
  ymax = 0.0257 # Test chamber radius
[]

[Variables]
  [./temp]
  [../]
[]

[Kernels]
  [./heat_conduction]
    type = HeatConduction
    variable = temp
  [../]
[]

[BCs]
  [./inlet_temperature]
    type = DirichletBC
    variable = temp
    boundary = left
    value = 350 # (K)
  [../]
  [./outlet_temperature]
    type = DirichletBC
    variable = temp
    boundary = right
    value = 300 # (K)
  [../]
[]

[Materials]
  [./steel]
    type = GenericConstantMaterial
    block = 0
    prop_names = thermal_conductivity
    prop_values = 18 # K: (W/m*K) from wikipedia @296K
  [../]
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  output_initial = true
  exodus = true
  print_pref_log = true
  print_linear_residuals = true
[]