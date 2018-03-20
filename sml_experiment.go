package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"
	"sync"

	"github.com/kzahedi/goent/discrete"
	pb "gopkg.in/cheggaaa/pb.v1"
)

var mutex sync.Mutex
var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")
var memprofile = flag.String("memprofile", "", "write mem profile to file")

var F_prime_xy_flag = flag.String("f_prime_xy", "diag", "Choice for the 'F_prime_xy' matrix for testing purposes.")

type McParameter struct {
	phi        float64
	psi        float64
	chi        float64
	mu         float64
	zeta       float64
	tau        float64
	mi         float64
	bins       int
	syci       int
	resolution int
}

type Indicator struct {
	a      int
	b      int
	c      int
	labmda float64
}

func lambda2(ind Indicator, aa, bb int) float64 {
	if ind.a == aa && ind.b == bb {
		return ind.labmda
	}
	return -ind.labmda //TODO: Why?
}

func lambda3(ind Indicator, aa, bb, cc int) float64 {
	if ind.a == aa && ind.b == bb && ind.c == cc {
		return ind.labmda
	}
	return -ind.labmda
}

func (ind *Indicator) String() string {
	return fmt.Sprintf("Indicator %d %d = %f", ind.a, ind.b, ind.labmda)
}

func check3DProbabilityDistribution(p [][][]float64) {
	sum := 0.0
	for x := 0; x < len(p); x++ {
		for y := 0; y < len(p[x]); y++ {
			for z := 0; z < len(p[x][y]); z++ {
				if math.IsNaN(p[x][y][z]) {
					panic("NaN")
				}
				sum += p[x][y][z]
			}
		}
	}
	if math.Abs(sum-1.0) > 0.0000001 {
		panic(fmt.Sprintf("Does not sum up to one %f", sum))
	}
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func vmap(v float64, bins int) float64 {
	return 2.0*float64(v)/float64(bins-1) - 1.0
}

func get_value2_pi(a, b float64, factor float64, bins int) float64 {
	S := 0.0
	var s float64
	for i := 0; i < bins; i += 1 {
		for j := 0; j < bins; j += 1{
			if *F_prime_xy_flag == "diag" {
				if i == j {
					s = 2*factor
				} else {
					s = 0.0
				}
			} else if *F_prime_xy_flag == "0000" {
				s = 0.0
			} else if *F_prime_xy_flag == "10101010" {
				s = 10.0
			} else if *F_prime_xy_flag == "flipped_diag" {
				if i == j {
					s = 0.0
				} else {
					s =  2.0*factor
				}
			} else {
				panic("unknown F_prime_xy_flag parameter.")
			}

			t1 := (2*float64(i) + 1.0) * a * math.Pi / float64(2 * bins) 
			t2 := (2*float64(j) + 1.0) * b * math.Pi / float64(2 * bins)
			S += s * math.Cos(t1) * math.Cos(t2)
		}
	}
	if a == 0 {
		S *= math.Sqrt(0.5)
	}
	if b == 0 {
		S *= math.Sqrt(0.5)
	}	
	return math.Sqrt(4 / float64(bins*bins)) * S

}

func get_value2(a, b float64, factor float64, bins int) float64 {
	return factor * vmap(a, bins) * vmap(b, bins)
}

func get_value3(a, b, c float64, factor float64, bins int) float64 {
	return factor * vmap(a, bins) * vmap(b, bins) * vmap(c, bins)
}

// func get_value3(a, b, c float64, factor float64, bins int) float64 {
// 	S := 0.0
// 	var s float64
// 	for i := 0; i < bins; i += 1 {
// 		for j := 0; j < bins; j += 1{
// 			for k := 0; k < bins; k += 1 {
// 				if i == j && j == k {
// 					s = 0.0 //* factor
// 				} else {
// 					s = factor
// 				}
// 				t1 := (2*float64(i) + 1.0) * a * math.Pi / float64(2 * bins) 
// 				t2 := (2*float64(j) + 1.0) * b * math.Pi / float64(2 * bins)
// 				t3 := (2*float64(k) + 1.0) * c * math.Pi / float64(2 * bins)
// 				S += s * math.Cos(t1) * math.Cos(t2) * math.Cos(t3)
// 			}
// 		}
// 	}
// 	if a == 0 {
// 		S *= math.Sqrt(0.5)
// 	}
// 	if b == 0 {
// 		S *= math.Sqrt(0.5)
// 	}	
// 	if c == 0 {
// 		S *= math.Sqrt(0.5)
// 	}
// 	return math.Sqrt(8 / float64(bins*bins*bins)) * S
// }

func f2s(v float64) string {
	return strconv.FormatFloat(v, 'f', -1, 64)
}

func getvalues(str string) []float64 {
	var r []float64
	comma_values := strings.Split(str, ",")

	if len(comma_values) == 1 {

		values := strings.Split(str, ":")
		start, err := strconv.ParseFloat(values[0], 64)
		check(err)

		end := start
		delta := start - end + 1.0

		if len(values) == 3 {
			delta, err = strconv.ParseFloat(values[1], 64)
			check(err)

			end, err = strconv.ParseFloat(values[2], 64)
			check(err)
		}

		for v := start; v <= end; v += delta {
			r = append(r, v)
		}
	} else {
		for _, v := range comma_values {
			f, _ := strconv.ParseFloat(v, 64)
			r = append(r, f)
		}
	}

	return r
}

func pw2_c_w1_a1(w2, w1, a1, bins int,
	phi, psi, chi float64,
	w2w1a1i, w2w1i, w2a1i []Indicator) float64 {
	z := 0.0
	n := 0.0
	for _, ind := range w2w1a1i {
		z += lambda3(ind, w2, w1, a1)
	}
	for _, ind := range w2w1i {
		z += lambda2(ind, w2, w1)
	}
	for _, ind := range w2a1i {
		z += lambda2(ind, w2, a1)
	}
	z = math.Exp(z)

	for w22 := 0; w22 < bins; w22++ {
		nn := 0.0

		for _, ind := range w2w1a1i {
			nn += lambda3(ind, w22, w1, a1)
		}
		for _, ind := range w2w1i {
			nn += lambda2(ind, w22, w1)
		}
		for _, ind := range w2a1i {
			nn += lambda2(ind, w22, a1)
		}
		n += math.Exp(nn)
	}

	return z / n
}

func pa1_c_s1(a1, s1, bins int, mu float64, a1s1i []Indicator) float64 {
	z := 0.0
	n := 0.0
	nn := 0.0
	for _, ind := range a1s1i {
		z += lambda2(ind, a1, s1)
	}
	z = math.Exp(z)

	for a11 := 0; a11 < bins; a11++ {
		nn = 0.0
		for _, ind := range a1s1i {
			nn += lambda2(ind, a11, s1)
		}
		n += math.Exp(nn)
	}
	return z / n
}

func ps1_c_w1(s1, w1, bins int, zeta float64, s1w1i []Indicator) float64 {
	z := 0.0
	n := 0.0
	nn := 0.0

	for _, ind := range s1w1i {
		z += lambda2(ind, s1, w1)
	}
	z = math.Exp(z)

	for s11 := 0; s11 < bins; s11++ {
		nn = 0.0
		for _, ind := range s1w1i {
			nn += lambda2(ind, s11, w1)
		}
		n += math.Exp(nn)
	}

	return z / n
}

func pw1(w1, bins int, tau float64) float64 {
	return 1.0 / float64(bins)
}

func generate_w2w1a1_indicators(chi float64, bins int) []Indicator {
	var r []Indicator
	for w2 := 0; w2 < bins; w2++ {
		for w1 := 0; w1 < bins; w1++ {
			for a1 := 0; a1 < bins; a1++ {
				r = append(r, Indicator{w2, w1, a1,
					get_value3(float64(w2), float64(w1), float64(a1), chi, bins)})
			}
		}
	}
	return r
}

func generate_w2w1_indicators(phi float64, bins int) []Indicator {
	var r []Indicator
	for w2 := 0; w2 < bins; w2++ {
		for w1 := 0; w1 < bins; w1++ {
			r = append(r, Indicator{w2, w1, -1,
				get_value2(float64(w2), float64(w1), phi, bins)})
		}
	}
	return r
}

func generate_w2a1_indicators(psi float64, bins int) []Indicator {
	var r []Indicator
	for w2 := 0; w2 < bins; w2++ {
		for a1 := 0; a1 < bins; a1++ {
			r = append(r, Indicator{w2, a1, -1,
				get_value2(float64(w2), float64(a1), psi, bins)})
		}
	}
	return r
}

func generate_a1s1_indicators(mu float64, bins int) []Indicator {
	var r []Indicator
	for a1 := 0; a1 < bins; a1++ {
		for s1 := 0; s1 < bins; s1++ {
			r = append(r, Indicator{a1, s1, -1,
				get_value2_pi(float64(a1), float64(s1), mu, bins)})
		}
	}
	return r
}

func generate_s1w1_indicators(zeta float64, bins int) []Indicator {
	var r []Indicator
	for s1 := 0; s1 < bins; s1++ {
		for w1 := 0; w1 < bins; w1++ {
			r = append(r, Indicator{s1, w1, -1,
				get_value2(float64(s1), float64(w1), zeta, bins)})
		}
	}
	return r
}

func calculate_MI_CA(p McParameter) McParameter {

	pw2w1 := discrete.Create2D(p.bins, p.bins)
	pw2a1 := discrete.Create2D(p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	v := 0.0

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					v = pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
						pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
						ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
						pw1(w1, p.bins, p.tau)
					pw2w1[w2][w1] += v
					pw2a1[w2][a1] += v
				}
			}
		}
	}

	r := discrete.MorphologicalComputationCW(pw2w1, pw2a1)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: r, bins: p.bins, tau: p.tau}
}

func calculate_CA(p McParameter) McParameter {

	ps2s1a1 := discrete.Create3D(p.bins, p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	v := 0.0

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					for s2 := 0; s2 < p.bins; s2++ {
						v = pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
							pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
							ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
							pw1(w1, p.bins, p.tau) *
							ps1_c_w1(s2, w2, p.bins, p.zeta, s1w1i) // this is for w' -> s'
						ps2s1a1[s2][s1][a1] += v
					}
				}
			}
		}
	}

	r := discrete.MorphologicalComputationIntrinsicCA(ps2s1a1, p.bins)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: r, bins: p.bins, tau: p.tau}
}

func calculate_CW(p McParameter) McParameter {

	ps2s1a1 := discrete.Create3D(p.bins, p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	v := 0.0

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					for s2 := 0; s2 < p.bins; s2++ {
						v = pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
							pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
							ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
							pw1(w1, p.bins, p.tau) *
							ps1_c_w1(s2, w2, p.bins, p.zeta, s1w1i) // this is for w' -> s'
						ps2s1a1[s2][s1][a1] += v
					}
				}
			}
		}
	}

	r := discrete.MorphologicalComputationIntrinsicCW(ps2s1a1)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: r, bins: p.bins, tau: p.tau}
}

func calculate_MI_WS(p McParameter) McParameter {

	pw2w1s1 := discrete.Create3D(p.bins, p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					pw2w1s1[w2][w1][s1] +=
						pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
							pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
							ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
							pw1(w1, p.bins, p.tau)
				}
			}
		}
	}

	r := discrete.MorphologicalComputationWS(pw2w1s1)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: r, bins: p.bins, tau: p.tau}
}

func calculate_MI_WA(p McParameter) McParameter {

	pw2w1a1 := discrete.Create3D(p.bins, p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					pw2w1a1[w2][w1][a1] +=
						pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
							pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
							ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
							pw1(w1, p.bins, p.tau)
				}
			}
		}
	}

	r := discrete.MorphologicalComputationWA(pw2w1a1)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: r, bins: p.bins, tau: p.tau}
}

func calculate_MI_W(p McParameter) McParameter {

	pw2w1a1 := discrete.Create3D(p.bins, p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					pw2w1a1[w2][w1][a1] +=
						pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
							pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
							ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
							pw1(w1, p.bins, p.tau)
				}
			}
		}
	}

	r := discrete.MorphologicalComputationW(pw2w1a1)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: r, bins: p.bins, tau: p.tau}
}

func calculate_MI_Wp(p McParameter) McParameter {

	pw2w1a1 := discrete.Create3D(p.bins, p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					pw2w1a1[w2][w1][a1] +=
						pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
							pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
							ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
							pw1(w1, p.bins, p.tau)
				}
			}
		}
	}

	r := discrete.MorphologicalComputationWp(pw2w1a1, p.syci, false)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: r, bins: p.bins, tau: p.tau}
}

func calculate_CI(p McParameter) McParameter {

	pw2w1a1 := discrete.Create3D(p.bins, p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					pw2w1a1[w2][w1][a1] +=
						pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
							pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
							ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
							pw1(w1, p.bins, p.tau)
				}
			}
		}
	}

	ci, _, _ := discrete.InformationDecomposition(pw2w1a1, p.resolution)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: ci, bins: p.bins, tau: p.tau}
}

func calculate_UI(p McParameter) McParameter {

	pw2w1a1 := discrete.Create3D(p.bins, p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					pw2w1a1[w2][w1][a1] +=
						pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
							pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
							ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
							pw1(w1, p.bins, p.tau)
				}
			}
		}
	}

	_, ui, _ := discrete.InformationDecomposition(pw2w1a1, p.resolution)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: ui, bins: p.bins, tau: p.tau}
}

func calculate_MI_MI(p McParameter) McParameter {

	pw2w1 := discrete.Create2D(p.bins, p.bins)
	pa1s1 := discrete.Create2D(p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	v := 0.0

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					v = pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
						pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
						ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
						pw1(w1, p.bins, p.tau)
					pw2w1[w2][w1] += v
					pa1s1[a1][s1] += v
				}
			}
		}
	}

	r := discrete.MorphologicalComputationMI(pw2w1, pa1s1)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: r, bins: p.bins, tau: p.tau}
}

func calculate_MI_A(p McParameter) McParameter {

	pw2a1w1 := discrete.Create3D(p.bins, p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					pw2a1w1[w2][a1][w1] +=
						pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
							pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
							ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
							pw1(w1, p.bins, p.tau)
				}
			}
		}
	}

	r := discrete.MorphologicalComputationA(pw2a1w1)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: r, bins: p.bins, tau: p.tau}
}

func calculate_MI_A_Prime(p McParameter) McParameter {

	pw2a1w1 := discrete.Create3D(p.bins, p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					pw2a1w1[w2][a1][w1] +=
						pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
							pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
							ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
							pw1(w1, p.bins, p.tau)
				}
			}
		}
	}

	r := 1.0 - discrete.MorphologicalComputationA(pw2a1w1)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: r, bins: p.bins, tau: p.tau}
}

func calculate_MI_SY(p McParameter) McParameter {

	pw2w1a1 := discrete.Create3D(p.bins, p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					pw2w1a1[w2][w1][a1] +=
						pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
							pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
							ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
							pw1(w1, p.bins, p.tau)
				}
			}
		}
	}

	r := discrete.MorphologicalComputationSY(pw2w1a1, p.syci, false)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: r, bins: p.bins, tau: p.tau}
}

func calculate_MI_SY_NID(p McParameter) McParameter {

	pw2w1a1 := discrete.Create3D(p.bins, p.bins, p.bins)
	w2w1a1i := generate_w2w1a1_indicators(p.chi, p.bins)
	w2w1i := generate_w2w1_indicators(p.phi, p.bins)
	w2a1i := generate_w2a1_indicators(p.psi, p.bins)
	a1s1i := generate_a1s1_indicators(p.mu, p.bins)
	s1w1i := generate_s1w1_indicators(p.zeta, p.bins)

	for w2 := 0; w2 < p.bins; w2++ {
		for w1 := 0; w1 < p.bins; w1++ {
			for a1 := 0; a1 < p.bins; a1++ {
				for s1 := 0; s1 < p.bins; s1++ {
					pw2w1a1[w2][w1][a1] +=
						pw2_c_w1_a1(w2, w1, a1, p.bins, p.phi, p.psi, p.chi, w2w1a1i, w2w1i, w2a1i) *
							pa1_c_s1(a1, s1, p.bins, p.mu, a1s1i) *
							ps1_c_w1(s1, w1, p.bins, p.zeta, s1w1i) *
							pw1(w1, p.bins, p.tau)
				}
			}
		}
	}

	r := discrete.MorphologicalComputationSyNid(pw2w1a1, p.syci)
	return McParameter{phi: p.phi, psi: p.psi, chi: p.chi, mu: p.mu, zeta: p.zeta, mi: r, bins: p.bins, tau: p.tau}
}

func main() {

	verbose := flag.Bool("v", false, "Print out the settings.")
	muStr := flag.String("mu", "0.0:0.5:1.0", "s -> a. can take list (1,2,3) or range with delta (0:0.1:1.0)")
	phiStr := flag.String("phi", "0.0:0.01:5", "w -> w'. Can take list (1,2,3) or range with delta (0:0.1:1.0)")
	psiStr := flag.String("psi", "0.0:0.01:5", "(a -> w'. Can take list (1,2,3) or range with delta (0:0.1:1.0)")
	chiStr := flag.String("chi", "0.0:1.25:5", "a,w -> w'. Can take list (1,2,3) or range with delta (0:0.1:1.0)")
	zetaStr := flag.String("zeta", "10.0", "w -> s. Can take list (1,2,3) or range with delta (0:0.1:1.0)")
	tauStr := flag.String("tau", "0", "p(w). Can take list (1,2,3) or range with delta (0:0.1:1.0)")
	mi := flag.String("mi", "MI_W", "available quantifications are: MI_W, MI_A, MI_A_Prime, MI_MI, MI_SY, MI_SY_NID, MI_CA, MI_WA, MI_WS, MI_Wp, CA, UI, CI")
	bins := flag.Int("b", 2, "Bins")
	syci := flag.Int("syci", 1000, "MI_SY convergence iterations")
	resolution := flag.Int("r", 100, "UI/CI resolution")
	cpus := flag.Int("cpus", 0, "Nr. of CPUs")
	output := flag.String("o", "./data", "output dir. Default: ./data")

	flag.Parse()

	if *verbose == true {
		fmt.Println(fmt.Sprintf("mu   (S -> A):    %s", *muStr))
		fmt.Println(fmt.Sprintf("phi  (W -> W'):   %s", *phiStr))
		fmt.Println(fmt.Sprintf("psi  (A -> W'):   %s", *psiStr))
		fmt.Println(fmt.Sprintf("chi  (A,W -> W'): %s", *chiStr))
		fmt.Println(fmt.Sprintf("tau  (W):         %s", *tauStr))
		fmt.Println(fmt.Sprintf("zeta (W -> S):    %s", *zetaStr))
		fmt.Println(fmt.Sprintf("Iterative scaling convergence iterations: %d", *syci))
		fmt.Println(fmt.Sprintf("UI/CI resolution: %d", *resolution))
	}

	if *cpuprofile != "" {
		fc, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(fc)
		defer pprof.StopCPUProfile()
	}

	if *cpus > 0 {
		runtime.GOMAXPROCS(*cpus)
	} else {
		runtime.GOMAXPROCS(runtime.NumCPU())
		*cpus = runtime.NumCPU()
	}

	mu := getvalues(*muStr)
	phi := getvalues(*phiStr)
	psi := getvalues(*psiStr)
	chi := getvalues(*chiStr)
	zeta := getvalues(*zetaStr)
	tau := getvalues(*tauStr)

	outputFiles := make(map[string]*os.File)

	fmt.Println(fmt.Sprintf("Creating directory %s", *output))
	rerr := os.Mkdir(*output, 0755)
	if rerr != nil {
		fmt.Println(fmt.Sprintf("Output directory %s already exists. Data will be added to that directory.", *output))
	}

	for _, _mu := range mu {
		for _, _chi := range chi {
			for _, _zeta := range zeta {
				for _, _tau := range tau {
					key := fmt.Sprintf("mu_%f_chi_%f_zeta_%f_tau_%f", _mu, _chi, _zeta, _tau)
					filename := fmt.Sprintf("%s/%s.csv", *output, key)
					f, err := os.Create(filename)
					check(err)
					outputFiles[key] = f
					defer f.Close()
					f.WriteString(fmt.Sprintf("# phi  (w -> w'):   %s\n", *phiStr))
					f.WriteString(fmt.Sprintf("# psi  (a -> w'):   %s\n", *psiStr))
					f.WriteString(fmt.Sprintf("# chi  (w,a -> w'): %f\n", _chi))
					f.WriteString(fmt.Sprintf("# mu   (s -> a):    %f\n", _mu))
					f.WriteString(fmt.Sprintf("# zeta (w -> s):    %f\n", _zeta))
					f.WriteString(fmt.Sprintf("# tau  (p(w):       %f (unused)\n", _tau))
					f.WriteString(fmt.Sprintf("# mi:               %s\n", *mi))
					f.WriteString(fmt.Sprintf("# bins:             %d\n", *bins))
					f.WriteString("#")
					if len(phi) > 1 {
						f.WriteString(" phi,")
					}
					if len(psi) > 1 {
						f.WriteString(" psi,")
					}
					if len(chi) > 1 {
						f.WriteString(" chi,")
					}
					// if len(mu) > 1 {
					// f.WriteString(" mu,")
					// }
					// if len(zeta) > 1 {
					// f.WriteString(" zeta,")
					// }
					// if len(tau) > 1 {
					// f.WriteString(" tau,")
					// }
					f.WriteString(fmt.Sprintf(" %s\n", *mi))
				}
			}
		}
	}

	var miFunc func(McParameter) McParameter

	switch *mi {
	case "MI_W":
		miFunc = calculate_MI_W
	case "MI_A":
		miFunc = calculate_MI_A
	case "MI_A_Prime":
		miFunc = calculate_MI_A_Prime
	case "MI_MI":
		miFunc = calculate_MI_MI
	case "MI_CA":
		miFunc = calculate_MI_CA
	case "CA":
		miFunc = calculate_CA
	case "CW":
		miFunc = calculate_CW
	case "MI_WS":
		miFunc = calculate_MI_WS
	case "MI_WA":
		miFunc = calculate_MI_WA
	case "MI_SY":
		miFunc = calculate_MI_SY
	case "MI_SY_NID":
		miFunc = calculate_MI_SY_NID
	case "MI_Wp":
		miFunc = calculate_MI_Wp
	case "CI":
		miFunc = calculate_CI
	case "UI":
		miFunc = calculate_UI
	default:
		panic(fmt.Sprintf("Unknown quantification given %s", *mi))
	}

	fmt.Println(fmt.Sprintf("Using %d cpus on %s with %d bins", *cpus, *mi, *bins))

	iterations := len(mu) * len(phi) * len(psi) * len(chi) * len(zeta) * len(tau)
	bar := pb.StartNew(iterations)

	send := make(chan McParameter, *cpus*2)
	ans := make(chan McParameter, *cpus*2)

	// start workers
	var wg sync.WaitGroup

	for i := 0; i < *cpus; i++ {
		wg.Add(1)
		go func(send <-chan McParameter, ans chan<- McParameter) {
			defer wg.Done()
			for p := range send {
				ans <- miFunc(p)
			}
		}(send, ans)
	}

	// start the jobs
	go func(send chan<- McParameter) {
		for _, vmu := range mu {
			for _, vphi := range phi {
				for _, vpsi := range psi {
					for _, vchi := range chi {
						for _, vzeta := range zeta {
							for _, vtau := range tau {
								send <- McParameter{phi: vphi, psi: vpsi, chi: vchi, mu: vmu, zeta: vzeta, tau: vtau, bins: *bins, mi: 0.0, syci: *syci, resolution: *resolution}
							}
						}
					}
				}
			}
		}
		close(send)
		wg.Wait()
		close(ans)
	}(send)

	s := ""
	for r := range ans {
		bar.Increment()
		s = ""
		if len(phi) > 1 {
			s = fmt.Sprintf("%s%.2f,", s, r.phi)
		}
		if len(psi) > 1 {
			s = fmt.Sprintf("%s%.2f,", s, r.psi)
		}
		// if len(chi) > 1 {
		// s = fmt.Sprintf("%s%.2f,", s, r.chi)
		// }
		// if len(mu) > 1 {
		// s = fmt.Sprintf("%s%.2f,", s, r.mu)
		// }
		// if len(zeta) > 1 {
		// s = fmt.Sprintf("%s%.2f,", s, r.zeta)
		// }
		if len(tau) > 1 {
			s = fmt.Sprintf("%s%.2f,", s, r.tau)
		}
		s = fmt.Sprintf("%s%.5f\n", s, r.mi)
		key := fmt.Sprintf("mu_%f_chi_%f_zeta_%f_tau_%f", r.mu, r.chi, r.zeta, r.tau)
		outputFiles[key].WriteString(s)
		// f.Sync()
	}

	bar.FinishPrint("Finished")
	if *memprofile != "" {
		fm, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.WriteHeapProfile(fm)
		fm.Close()
		return
	}
}
