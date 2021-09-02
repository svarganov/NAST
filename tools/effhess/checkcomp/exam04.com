&programs
Molpro = .true.
&end

&molpro_input_files
low_spin_grad = 'HOHON_anion_singlet_grad_molpro.out'
high_spin_grad = 'HOHON_anion_triplet_grad_molpro.out'
low_spin_hess = 'HOHON_anion_singlet_hess_molpro.out'
high_spin_hess = 'HOHON_anion_triplet_hess_molpro.out'
&end
