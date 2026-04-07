&programs
ORCA = .true.
&end

&orca_input_files
low_spin_grad_orca = 'cyclopropene_S.out'
high_spin_grad_orca = 'cyclopropene_T.out'
low_spin_hess_orca = 'cyclopropene_S.hess'
high_spin_hess_orca = 'cyclopropene_T.hess'
&end
