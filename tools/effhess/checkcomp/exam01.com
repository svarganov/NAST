&programs
GAMESS = .true.
Molpro = .false.
&end

&gamess_input_files
low_spin_dat = 'methoxy_cation_singlet_grad_and_hess.dat'
high_spin_dat = 'methoxy_cation_triplet_grad_and_hess.dat'
input_geometry = 'methoxy_cation_singlet_grad_and_hess.inp'
&end
