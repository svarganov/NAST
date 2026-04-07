&programs
ORCA = .true.
&end

&orca_input_file
hess_orca = 'diketone_S1.hess'
&end

&geometries
initial_geometry = 'diketone_S1.xyz'
final_geometry = 'diketone_T3.xyz'
&end
