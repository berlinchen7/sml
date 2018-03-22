#!/bin/bash

# go run sml_dct.go -zeta 10.0 -psi 0.0:0.1:5.01 -mu 0.0:0.3:0.91 -phi 0.0:0.1:5.01 -chi 0.0:1.0:10.01 -o ./data_binary 
# go run sml_dct.go -zeta 10.0 -psi 0.0:0.1:5.01 -mu 0.0:0.3:0.91 -phi 0.0:0.1:5.01 -chi 0.0:1.0:10.01 -mi MI_SY -o ./data_binary_MSY 

go run sml_dct.go -zeta 10.0 -psi 0.0:0.1:5.01 -mu 0.0:0.3:0.91 -phi 0.0:0.1:5.01 -chi 0.0:1.0:10.01 -F_prime_config not_inv -o ./data_binary_none_inv 
go run sml_dct.go -zeta 10.0 -psi 0.0:0.1:5.01 -mu 0.0:0.3:0.91 -phi 0.0:0.1:5.01 -chi 0.0:1.0:10.01 -mi MI_SY -F_prime_config not_inv -o ./data_binary_MSY_none_inv 


echo done!