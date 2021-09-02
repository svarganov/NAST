&programs
GAMESS = .true.
Molpro = .false.
&end

&gamess_input_files
low_spin_dat = 'cyclopropene_singlet_grad_and_hess.dat'
high_spin_dat = 'cyclopropene_triplet_grad_and_hess.dat'
input_geometry = 'cyclopropene_triplet_grad_and_hess.inp'
&end

