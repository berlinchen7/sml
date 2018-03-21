package main

import "fmt"

func main() {

	const bins = 3

	var s [bins][bins][bins]float64

	for k := 0; k < bins; k += 1 {
		for i := 0; i < bins; i += 1 {
			for j := 0; j < bins; j += 1 {
				if i + j == k {
					s[i][j][k] = 1.0
				} else {
					s[i][j][k] = 0.0
				}
			}
		}
	}

	fmt.Println(s)

}