CDF      
      
len_string     !   len_line   Q   four      	time_step          len_name   !   num_dim       	num_nodes         num_elem      
num_el_blk        num_node_sets         num_el_in_blk1        num_nod_per_el1       num_el_in_blk2        num_nod_per_el2       num_nod_ns1       num_nod_ns2       num_nod_ns3       num_nod_var       num_elem_var      num_info  �         api_version       @�
=   version       @�
=   floating_point_word_size            	file_size               int64_status             title         solid_mech_truss_out.e     maximum_name_length                    
time_whole                            �    	eb_status                             t   eb_prop1               name      ID              |   	ns_status         	                    �   ns_prop1      	         name      ID              �   coordx                            �   coordy                            �   coordz                            �   eb_names                       D      �   ns_names      	                 d      	(   
coor_names                         d      	�   connect1      
            	elem_type         EDGE2               	�   elem_num_map                          	�   connect2                  	elem_type         EDGE2               
    node_ns1                          
   node_ns2                          
   node_ns3                          
   vals_nod_var1                                �   vals_nod_var2                                �    vals_nod_var3                                �8   name_nod_var                       d      
   name_elem_var                          D      
x   vals_elem_var1eb1            
                    �P   vals_elem_var2eb1            
                    �X   vals_elem_var1eb2                                �`   vals_elem_var2eb2                                �h   elem_var_tab                             
�   info_records                      4      
�                                      ?�      @                                                                                                                                                                                                                                                                                                                                                              disp_x                           disp_y                           disp_z                            axial_stress                     e_over_l                                       ####################                                                             # Created by MOOSE #                                                             ####################                                                             ### Command Line Arguments ###                                                   -i                                                                               solid_mech_truss.i                                                                                                                                                ### Version Info ###                                                             Framework Information:                                                           SVN Revision:            15a3143506166ec5cc82f515f65c1b31cc10ab57                PETSc Version:           3.3.0                                                   Current Time:            Mon Sep 23 10:14:33 2013                                Executable Timestamp:    Mon Sep 23 10:14:13 2013                                                                                                                                                                                                  ### Input File ###                                                                                                                                                []                                                                                 block                          =                                                 coord_type                     = XYZ                                             fe_cache                       = 0                                               name                           = 'MOOSE Problem'                                 type                           = FEProblem                                       order                          = AUTO                                            active_bcs                     =                                                 active_kernels                 =                                                 inactive_bcs                   =                                                 inactive_kernels               =                                                 start                          = 0                                               decomposition                  =                                                 dimNearNullSpace               = 0                                               dimNullSpace                   = 0                                               long_name                      =                                                 solve                          = 1                                             []                                                                                                                                                                [AuxKernels]                                                                                                                                                        [./axial_stress]                                                                   type                         = MaterialRealAux                                   block                        =                                                   boundary                     =                                                   execute_on                   = residual                                          factor                       = 1                                                 long_name                    = AuxKernels/axial_stress                           offset                       = 0                                                 property                     = axial_stress                                      use_displaced_mesh           = 0                                                 variable                     = axial_stress                                    [../]                                                                                                                                                             [./e_over_l]                                                                       type                         = MaterialRealAux                                   block                        =                                                   boundary                     =                                                   execute_on                   = residual                                          factor                       = 1                                                 long_name                    = AuxKernels/e_over_l                               offset                       = 0                                                 property                     = e_over_l                                          use_displaced_mesh           = 0                                                 variable                     = e_over_l                                        [../]                                                                          []                                                                                                                                                                [AuxVariables]                                                                                                                                                      [./axial_stress]                                                                   family                       = MONOMIAL                                          initial_condition            = 0                                                 order                        = CONSTANT                                          initial_from_file_timestep   = 2                                                 initial_from_file_var        =                                                 [../]                                                                                                                                                             [./e_over_l]                                                                       family                       = MONOMIAL                                          initial_condition            = 0                                                 order                        = CONSTANT                                          initial_from_file_timestep   = 2                                                 initial_from_file_var        =                                                 [../]                                                                          []                                                                                                                                                                [BCs]                                                                                                                                                               [./fixx1]                                                                          boundary                     = 1                                                 type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = disp_x                                            long_name                    = BCs/fixx1                                         value                        = 0                                               [../]                                                                                                                                                             [./fixx2]                                                                          boundary                     = 2                                                 type                         = FunctionDirichletBC                               use_displaced_mesh           = 0                                                 variable                     = disp_x                                            function                     = x2                                                long_name                    = BCs/fixx2                                       [../]                                                                                                                                                             [./fixx3]                                                                          boundary                     = 3                                                 type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = disp_x                                            long_name                    = BCs/fixx3                                         value                        = 0                                               [../]                                                                                                                                                             [./fixy1]                                                                          boundary                     = 1                                                 type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = disp_y                                            long_name                    = BCs/fixy1                                         value                        = 0                                               [../]                                                                                                                                                             [./fixy2]                                                                          boundary                     = 2                                                 type                         = FunctionDirichletBC                               use_displaced_mesh           = 0                                                 variable                     = disp_y                                            function                     = y2                                                long_name                    = BCs/fixy2                                       [../]                                                                                                                                                             [./fixy3]                                                                          boundary                     = 3                                                 type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = disp_y                                            long_name                    = BCs/fixy3                                         value                        = 0                                               [../]                                                                                                                                                             [./fixz1]                                                                          boundary                     = 1                                                 type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = disp_z                                            long_name                    = BCs/fixz1                                         value                        = 0                                               [../]                                                                                                                                                             [./fixz2]                                                                          boundary                     = 2                                                 type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = disp_z                                            long_name                    = BCs/fixz2                                         value                        = 0                                               [../]                                                                                                                                                             [./fixz3]                                                                          boundary                     = 3                                                 type                         = DirichletBC                                       use_displaced_mesh           = 0                                                 variable                     = disp_z                                            long_name                    = BCs/fixz3                                         value                        = 0                                               [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      l_abs_step_tol                 = -1                                              l_max_its                      = 10000                                           l_tol                          = 1e-05                                           line_search                    = default                                         nl_abs_step_tol                = 1e-50                                           nl_abs_tol                     = 1e-10                                           nl_max_funcs                   = 10000                                           nl_max_its                     = 15                                              nl_rel_step_tol                = 1e-50                                           nl_rel_tol                     = 1e-06                                           no_fe_reinit                   = 0                                               petsc_options                  = '-snes_mf_operator -ksp_monitor'                petsc_options_iname            = '-pc_type -snes_type -snes_ls -snes_linese... arch_type -ksp_gmres_restart'                                                      petsc_options_value            = 'jacobi ls basic basic 101'                     solve_type                     =                                                 type                           = Transient                                       _fe_problem                    = 0x7ffd4209dc00                                  abort_on_solve_fail            = 0                                               dt                             = 1                                               dtmax                          = 1e+30                                           dtmin                          = 2e-14                                           end_time                       = 3                                               long_name                      = Executioner                                     n_startup_steps                = 0                                               num_steps                      = 3                                               predictor_scale                =                                                 restart_file_base              =                                                 scheme                         = implicit-euler                                  ss_check_tol                   = 1e-08                                           ss_tmin                        = 0                                               start_time                     = 0                                               sync_times                     = -1.79769e+308                                   time_period_ends               =                                                 time_period_starts             =                                                 time_periods                   =                                                 timestep_tolerance             = 2e-14                                           trans_ss_check                 = 0                                               use_multiapp_dt                = 0                                             []                                                                                                                                                                [Functions]                                                                                                                                                         [./x2]                                                                             type                         = PiecewiseLinear                                   axis                         =                                                   long_name                    = Functions/x2                                      scale_factor                 = 1                                                 x                            = '0 1 2 3'                                         y                            = '0 0.5 1 1'                                     [../]                                                                                                                                                             [./y2]                                                                             type                         = PiecewiseLinear                                   axis                         =                                                   long_name                    = Functions/y2                                      scale_factor                 = 1                                                 x                            = '0 1 2 3'                                         y                            = '0 0 0.5 1'                                     [../]                                                                          []                                                                                                                                                                [Materials]                                                                                                                                                         [./density]                                                                        type                         = Density                                           block                        = '1 2'                                             boundary                     =                                                   density                      = 1                                                 disp_r                       =                                                   disp_x                       =                                                   disp_y                       =                                                   disp_z                       =                                                   long_name                    = Materials/density                               [../]                                                                                                                                                             [./dummy]                                                                          type                         = Elastic                                           active_crack_planes          =                                                   appended_property_name       =                                                   block                        = '1 2'                                             boundary                     =                                                   bulk_modulus                 =                                                   cracking_release             = abrupt                                            cracking_residual_stress     = 0                                                 cracking_stress              = 0                                                 disp_r                       =                                                   disp_x                       = disp_x                                            disp_y                       = disp_y                                            disp_z                       = disp_z                                            formulation                  =                                                   increment_calculation        = RashidApprox                                      lambda                       =                                                   large_strain                 = 0                                                 long_name                    = Materials/dummy                                   max_cracks                   = 3                                                 poissons_ratio               = 0.2                                               poissons_ratio_function      =                                                   shear_modulus                =                                                   stress_free_temperature      =                                                   temp                         =                                                   thermal_expansion            =                                                   thermal_expansion_function   =                                                   youngs_modulus               = 1e+06                                             youngs_modulus_function      =                                                 [../]                                                                                                                                                             [./linelast]                                                                       type                         = TrussMaterial                                     block                        = '1 2'                                             boundary                     =                                                   disp_x                       = disp_x                                            disp_y                       = disp_y                                            disp_z                       = disp_z                                            long_name                    = Materials/linelast                                t_ref                        = 0                                                 thermal_expansion            = 0                                                 youngs_modulus               = 1e+06                                           [../]                                                                          []                                                                                                                                                                [Mesh]                                                                             displacements                  =                                                 block_id                       =                                                 block_name                     =                                                 boundary_id                    =                                                 boundary_name                  =                                                 centroid_partitioner_direction =                                                 construct_side_list_from_node_list = 0                                           ghosted_boundaries             =                                                 ghosted_boundaries_inflation   =                                                 partitioner                    = metis                                           patch_size                     = 40                                              second_order                   = 0                                               type                           = FileMesh                                        uniform_refine                 = 0                                               dim                            = 3                                               distribution                   = DEFAULT                                         file                           = twoTrussLine.e                                  long_name                      = Mesh                                            nemesis                        = 0                                               skip_partitioning              = 0                                             []                                                                                                                                                                [Output]                                                                           elemental_as_nodal             = 0                                               exodus                         = 1                                               exodus_inputfile_output        = 1                                               file_base                      = solid_mech_truss_out                            gmv                            = 0                                               gnuplot_format                 = ps                                              hidden_variables               =                                                 interval                       = 1                                               iteration_plot_start_time      = 1.79769e+308                                    max_pps_rows_screen            = 15                                              nemesis                        = 0                                               num_restart_files              = 0                                               output_displaced               = 0                                               output_es_info                 = 1                                               output_if_base_contains        =                                                 output_initial                 = 1                                               output_solution_history        = 0                                               output_variables               =                                                 perf_log                       = 1                                               position                       =                                                 postprocessor_csv              = 0                                               postprocessor_gnuplot          = 0                                               postprocessor_screen           = 1                                               pps_fit_to_screen              = ENVIRONMENT                                     print_linear_residuals         = 0                                               screen_interval                = 1                                               show_setup_log_early           = 0                                               tecplot                        = 0                                               tecplot_binary                 = 0                                               time_interval                  =                                                 vtk                            = 0                                               xda                            = 0                                             []                                                                                                                                                                [SolidMechanics]                                                                                                                                                    [./solid]                                                                          appended_property_name       =                                                   block                        =                                                   diag_save_in_disp_r          =                                                   diag_save_in_disp_x          =                                                   diag_save_in_disp_y          =                                                   diag_save_in_disp_z          =                                                   disp_r                       =                                                   disp_x                       = disp_x                                            disp_y                       = disp_y                                            disp_z                       = disp_z                                            save_in_disp_r               =                                                   save_in_disp_x               =                                                   save_in_disp_y               =                                                   save_in_disp_z               =                                                   temp                         =                                                   use_displaced_mesh           = 1                                               [../]                                                                          []                                                                                                                                                                [Variables]                                                                                                                                                         [./disp_x]                                                                         block                        =                                                   family                       = LAGRANGE                                          initial_condition            = 0                                                 order                        = FIRST                                             scaling                      = 1                                                 initial_from_file_timestep   = 2                                                 initial_from_file_var        =                                                 [../]                                                                                                                                                             [./disp_y]                                                                         block                        =                                                   family                       = LAGRANGE                                          initial_condition            = 0                                                 order                        = FIRST                                             scaling                      = 1                                                 initial_from_file_timestep   = 2                                                 initial_from_file_var        =                                                 [../]                                                                                                                                                             [./disp_z]                                                                         block                        =                                                   family                       = LAGRANGE                                          initial_condition            = 0                                                 order                        = FIRST                                             scaling                      = 1                                                 initial_from_file_timestep   = 2                                                 initial_from_file_var        =                                                 [../]                                                                          []                                                                                                                                                                                                 ?�              ?�                                                             A��    A.��    ���   A.��    @               ?�   }h                ?�   2��                                A02��1��A.��    ����ÍA.��    @              ?�                      ?�   )��                                A2�c�O:�A.��    ?C��   A.��    