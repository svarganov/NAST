&programs
ORCA = .true.
&end

&orca_input_file
hess_orca = 'TTM_1Cz_An_2D1S0.hess'
&end

&geometries
initial_geometry = 'TTM_1Cz_An_2D0T1.xyz'
final_geometry = 'TTM_1Cz_An_2D1S0.xyz'
&end
