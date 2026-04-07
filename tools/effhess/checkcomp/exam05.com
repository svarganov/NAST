&programs
QChem = .true.
&end

&qchem_input_files
geometry_in = 'cyclopropene_mecp_hess_S.in'
low_spin_grad_fchk = 'cyclopropene_mecp_force_S.fchk'
high_spin_grad_fchk = 'cyclopropene_mecp_force_T.fchk'
low_spin_hess_fchk = 'cyclopropene_mecp_hess_S.fchk'
high_spin_hess_fchk = 'cyclopropene_mecp_hess_T.fchk'
&end
