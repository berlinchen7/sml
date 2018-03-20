#!/bin/bash
# Create MI_W output based on various choices of T_{u, v}.

go run sml.go -zeta 10.0 -psi 0.0 -mu 0.0:0.25:1.01 -phi 5.0 -chi 0.0 -o ./data_orig
go run sml_experiment.go -zeta 10.0 -psi 0.0 -mu 0.0:0.25:1.01 -phi 5.0 -chi 0.0 -f_prime_xy diag -o ./data_diag 
go run sml_experiment.go -zeta 10.0 -psi 0.0 -mu 0.0:0.25:1.01 -phi 5.0 -chi 0.0 -f_prime_xy 0000 -o ./data_0000
go run sml_experiment.go -zeta 10.0 -psi 0.0 -mu 0.0:0.25:1.01 -phi 5.0 -chi 0.0 -f_prime_xy 10101010 -o ./data_10101010
go run sml_experiment.go -zeta 10.0 -psi 0.0 -mu 0.0:0.25:1.01 -phi 5.0 -chi 0.0 -f_prime_xy flipped_diag -o ./data_flipped_diag

echo done!