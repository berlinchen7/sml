package main

import "fmt"

func main() {

	const bins = 3

	var s [bins][bins][bins]float64

	for k := 0; k < bins; k += 1 {
		for i := 0; i < bins; i += 1 {
			for j := 0; j < bins; j += 1 {
				if ((bins - 1) - i + j == k - 1) || ((bins - 1) - i + j == k + (bins - 1) ) {
					s[i][j][k] = 1.0
				} else {
					s[i][j][k] = 0.0
				}
			}
		}
	}

	for a := 0; a < bins; a += 1 {
		for b := 0; b < bins; b += 1{
			fmt.Println(s[a][b][2])
		}
	}

}