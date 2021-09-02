&programs
GAMESS = .true.
&end

&gamess_input_files
low_spin_dat = 'Ni_dpp_singlet_grad_and_hess.dat'
high_spin_dat = 'Ni_dpp_triplet_grad_and_hess.dat'
input_geometry = 'Ni_dpp_singlet_grad_and_hess.inp'
&end
