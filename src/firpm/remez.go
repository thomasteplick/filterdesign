/**************************************************************************
 * remez.go
 *
 * Parks-McClellan Algorithm for FIR filter design.  The particular criterion
 * used in this design procedure is the so-called minimax or Chebyshev criterion
 * where within the frequency intervals of interest (the passband and stopband
 * for a lowpass filter) we seek a frequency response that minimizes the maximum
 * weighted approximation error.  The Parks-McClellan algorithm for optimum
 * equiripple approximation of FIR filters can be used to design a wide variety
 * of FIR filters, such as lowpass, highpass, bandpass, bandstop, differentiator,
 * and hilbert transform.

 * Reference Discrete-Time Signal Processing, A. Oppenheim & R. Schafer, 1989
 * This code was translated from github.com/janovetz/remez-exchange C code.
 * See the copyright statement there.  It is GNU Library General Public License
 * as published by the the Free Software Foundation.
 *************************************************************************/

package firpm

import (
	"fmt"
	"math"
)

const (
	pi             float64 = math.Pi
	twoPi          float64 = 2.0 * pi
	maxIterations          = 40
	minGridDensity         = 16
)

// filter symmetry
const (
	negative = iota // 0
	positive        // 1
)

// filter types
const (
	lowpass        = iota // 0, not used
	bandpass              // 1
	differentiator        // 2
	hilbert               // 3
	highpass              // 4, not used
)

// Parks-McClellan FIR state variables
type Firpm struct {
	Grid      []float64
	D         []float64
	W         []float64
	E         []float64
	Extremals []int
	taps      []float64
	x         []float64
	y         []float64
	ad        []float64
}

// newFirpm returns a pointer to a Firpm struct
func newFirpm(gridsize int, r int) *Firpm {

	return &Firpm{Grid: make([]float64, gridsize),
		D:         make([]float64, gridsize),
		W:         make([]float64, gridsize),
		E:         make([]float64, gridsize),
		Extremals: make([]int, r+1),
		taps:      make([]float64, r+1),
		x:         make([]float64, r+1),
		y:         make([]float64, r+1),
		ad:        make([]float64, r+1),
	}
}

// createDenseGrid creates the dense grid of frequencies from the specified bands.
func (firpm *Firpm) createDenseGrid(r int, numtaps int, gridsize int, griddensity int,
	bands []float64, des []float64, weight []float64, sym int) {
	var (
		lowf  float64
		highf float64
		grid0 float64
	)

	delf := 0.5 / (float64(griddensity * r))
	numband := len(bands) / 2

	/*
	* For differentiator, hilbert,
	*   symmetry is odd and Grid[0] = max(delf, bands[0])
	 */
	if sym == negative && delf > bands[0] {
		grid0 = delf
	} else {
		grid0 = bands[0]
	}

	j := 0
	for band := 0; band < numband; band++ {
		if band == 0 {
			lowf = grid0
		} else {
			lowf = bands[2*band]
		}
		highf = bands[2*band+1]
		k := (int)((highf-lowf)/delf + 0.5) /* .5 for rounding */
		for i := 0; i < k; i++ {
			firpm.D[j] = des[2*band] + float64(i)*(des[2*band+1]-des[2*band])/float64((k-1))
			firpm.W[j] = weight[band]
			firpm.Grid[j] = lowf
			lowf += delf
			j++
		}
		firpm.Grid[j-1] = highf
	}

	/*
	* Similar to above, if odd symmetry, last grid point can't be .5
	*  - but, if there are even taps, leave the last grid point at .5
	 */
	if (sym == negative) &&
		(firpm.Grid[gridsize-1] > (0.5 - delf)) &&
		numtaps%2 != 0 {
		firpm.Grid[gridsize-1] = 0.5 - delf
	}
}

// initialGuess places extremal frequencies evenly throughout the dense grid
func (firpm *Firpm) initialGuess(r int, gridsize int) {
	for i := 0; i <= r; i++ {
		firpm.Extremals[i] = i * (gridsize - 1) / r
	}
}

// calcParams performs Oppenheim & Schafer 7.131, 7.132, and 7.133b to get (7.133a)
func (firpm *Firpm) calcParams(r int) {
	// delta (7.131), b (7.132) is ad[], C (7.133b) is y[]
	// These are used to get A(e^jomega) (7.133a) and also E(omega) (7.113)

	/*
	 * Find x[]
	 */
	for i := 0; i <= r; i++ {
		firpm.x[i] = math.Cos(twoPi * firpm.Grid[firpm.Extremals[i]])
	}
	/*
	 * Calculate ad[]  - Oppenheim & Schafer eq 7.132
	 */
	ld := (r-1)/15 + 1 /* Skips around to avoid round errors */
	for i := 0; i <= r; i++ {
		denom := 1.0
		xi := firpm.x[i]
		for j := 0; j < ld; j++ {
			for k := j; k <= r; k += ld {
				if k != i {
					denom *= 2.0 * (xi - firpm.x[k])
				}
			}
		}
		if math.Abs(denom) < 0.00001 {
			denom = 0.00001
		}
		firpm.ad[i] = 1.0 / denom
	}

	/*
	 * Calculate delta  - Oppenheim & Schafer eq 7.131
	 */
	numer := 0.0
	denom := 0.0
	sign := 1.0
	for i := 0; i <= r; i++ {
		numer += firpm.ad[i] * firpm.D[firpm.Extremals[i]]
		denom += sign * firpm.ad[i] / firpm.W[firpm.Extremals[i]]
		sign = -sign
	}
	delta := numer / denom
	sign = 1.0

	/*
	 * Calculate y[]  - Oppenheim & Schafer eq 7.133b
	 */
	for i := 0; i <= r; i++ {
		firpm.y[i] = firpm.D[firpm.Extremals[i]] - sign*delta/firpm.W[firpm.Extremals[i]]
		sign = -sign
	}
}

// computeA uses calcParams to perform Oppenheim & Schafer (7.133a) A(e^jomega)
func (firpm *Firpm) computeA(freq float64, r int) float64 {
	// Calculates the actual filter response at a given frequency A(e^jomega)
	denom := 0.0
	numer := 0.0
	xc := math.Cos(twoPi * freq)
	for i := 0; i <= r; i++ {
		c := xc - firpm.x[i]
		if math.Abs(c) < 1.0e-7 {
			numer = firpm.y[i]
			denom = 1
			break
		}
		c = firpm.ad[i] / c
		denom += c
		numer += c * firpm.y[i]
	}
	return numer / denom
}

// calcError calculates the Error function from the desired frequency response (7.113)
func (firpm *Firpm) calcError(r int, gridsize int) {
	// E(omega) Oppenheim & Schafer (7.113)

	for i := 0; i < gridsize; i++ {
		A := firpm.computeA(firpm.Grid[i], r)
		firpm.E[i] = firpm.W[i] * (firpm.D[i] - A)
	}
}

// search finds the maxima/minima of the error curve.
func (firpm *Firpm) search(r int, gridsize int) error {

	foundExt := make([]int, 2*r)
	k := 0

	/*
	* Check for extremum at 0.
	 */
	if ((firpm.E[0] > 0.0) && (firpm.E[0] > firpm.E[1])) ||
		((firpm.E[0] < 0.0) && (firpm.E[0] < firpm.E[1])) {
		foundExt[k] = 0
		k++
	}
	/*
	* Check for extrema inside dense grid
	 */
	for i := 1; i < gridsize-1; i++ {
		if ((firpm.E[i] >= firpm.E[i-1]) && (firpm.E[i] > firpm.E[i+1]) && (firpm.E[i] > 0.0)) ||
			((firpm.E[i] <= firpm.E[i-1]) && (firpm.E[i] < firpm.E[i+1]) && (firpm.E[i] < 0.0)) {
			// PAK: we sometimes get too many extremal frequencies
			if k >= 2*r {
				return fmt.Errorf("Remez: too many extremals->cannot continue")
			}
			foundExt[k] = i
			k++
		}
	}

	/*
	* Check for extremum at 0.5
	 */
	j := gridsize - 1
	if ((firpm.E[j] > 0.0) && (firpm.E[j] > firpm.E[j-1])) ||
		((firpm.E[j] < 0.0) && (firpm.E[j] < firpm.E[j-1])) {
		if k >= 2*r {
			return fmt.Errorf("Remez: too many extremals->cannot continue")
		}
		foundExt[k] = j
		k++
	}

	// PAK: we sometimes do not get enough extremal frequencies
	if k < r+1 {
		return fmt.Errorf("Remez:  insufficient extremals->cannot continue")
	}

	/*
	* Remove extra extremals
	 */
	extra := k - (r + 1)

	var (
		up  int
		alt int
	)
	for extra > 0 {
		if firpm.E[foundExt[0]] > 0.0 {
			up = 1 /* first one is a maxima */
		} else {
			up = 0 /* first one is a minima */
		}
		l := 0
		alt = 1
		for j := 1; j < k; j++ {
			if math.Abs(firpm.E[foundExt[j]]) < math.Abs(firpm.E[foundExt[l]]) {
				l = j /* new smallest error. */
			}
			if up == 1 && (firpm.E[foundExt[j]] < 0.0) {
				up = 0 /* switch to a minima */
			} else if (up == 0) && (firpm.E[foundExt[j]] > 0.0) {
				up = 1 /* switch to a maxima */
			} else {
				alt = 0
				// PAK: break now and you will delete the smallest overall
				// extremal.  If you want to delete the smallest of the
				// pair of non-alternating extremals, then you must do:
				//
				/*
					if math.Abs(firpm.E[foundExt[j]]) < math.Abs(firpm.E[foundExt[j-1]]) {
						l = j
					} else {
						l = j - 1
					}
				*/
				break /* Ooops, found two non-alternating */
			} /* extrema.  Delete smallest of them */
		} /* if the loop finishes, all extrema are alternating */

		/*
		* If there's only one extremal and all are alternating,
		* delete the smallest of the first/last extremals.
		 */
		if alt == 1 && (extra == 1) {
			if math.Abs(firpm.E[foundExt[k-1]]) < math.Abs(firpm.E[foundExt[0]]) {
				/* Delete last extremal */
				l = k - 1
				// PAK: changed from l = foundExt[k-1];
			} else {
				/* Delete first extremal */
				l = 0
				// PAK: changed from l = foundExt[0];
			}
		}
		for j = l; j < k-1; j++ { /* Loop that does the deletion */
			foundExt[j] = foundExt[j+1]
			//  assert(foundExt[j]<gridsize);
		}
		k--
		extra--
	}

	for i := 0; i <= r; i++ {
		// assert(foundExt[i]<gridsize);
		firpm.Extremals[i] = foundExt[i] /* Copy found extremals to Extremals[] */
	}

	return nil
}

// freqSample is a frequency sampling algorithm to determine the impulse response h[]
// from A's found in computeA
func (firpm *Firpm) freqSample(N int, A []float64, h []float64, symm int) {

	var M float64 = (float64(N) - 1.0) / 2.0
	if symm == positive {
		if N%2 != 0 {
			for n := 0; n < N; n++ {
				val := A[0]
				x := twoPi * (float64(n) - M) / float64(N)
				for k := 1; k <= int(M); k++ {
					val += 2.0 * A[k] * math.Cos(x*float64(k))
				}
				h[n] = val / float64(N)
			}
		} else {
			for n := 0; n < N; n++ {
				val := A[0]
				x := twoPi * (float64(n) - M) / float64(N)
				for k := 1; k <= (N/2 - 1); k++ {
					val += 2.0 * A[k] * math.Cos(x*float64(k))
				}
				h[n] = val / float64(N)
			}
		}
	} else {
		if N%2 != 0 {
			for n := 0; n < N; n++ {
				val := 0.0
				x := twoPi * (float64(n) - M) / float64(N)
				for k := 1; k <= int(M); k++ {
					val += 2.0 * A[k] * math.Sin(x*float64(k))
				}
				h[n] = val / float64(N)
			}
		} else {
			for n := 0; n < N; n++ {
				val := A[N/2] * math.Sin(pi*(float64(n)-M))
				x := twoPi * (float64(n) - M) / float64(N)
				for k := 1; k <= (N/2 - 1); k++ {
					val += 2.0 * A[k] * math.Sin(x*float64(k))
				}
				h[n] = val / float64(N)
			}
		}
	}
}

// isDone checks to see if the error function is small enough to consider the result
// to have converged
func (firpm *Firpm) isDone(r int) bool {

	min := math.Abs(firpm.E[firpm.Extremals[0]])
	max := math.Abs(firpm.E[firpm.Extremals[0]])
	for i := 1; i <= r; i++ {
		current := math.Abs(firpm.E[firpm.Extremals[i]])
		if current < min {
			min = current
		}
		if current > max {
			max = current
		}
	}
	return (((max - min) / max) < 0.0001)
}

// Remez calculates the optimal (in the Chebyshev/minimax sense) FIR filter using the Parks-McClellan algorithm
func Remez(order int, f []float64, a []float64, w []float64,
	ftype int, density int) ([]float64, error) {

	var (
		symmetry int
		firpm    *Firpm
		h        []float64
		numband  int = len(f) / 2
	)

	if ftype == bandpass {
		symmetry = positive
	} else {
		symmetry = negative
	}

	// number of coefficients in impulse response of FIR filter
	numtaps := order + 1
	r := numtaps / 2 /* number of extrema */
	if (numtaps%2) != 0 && (symmetry == positive) {
		r++
	}

	// the FIR filter coefficients to be returned
	h = make([]float64, numtaps)
	h[0] = 32

	/*
	* Predict dense grid size in advance for memory allocation,
	* .5 is used so we round up, not truncate
	 */
	gridsize := 0
	for i := 0; i < numband; i++ {
		gridsize += (int)(2.0*float64(r)*float64(density)*(f[2*i+1]-f[2*i]) + .5)
	}
	if symmetry == negative {
		gridsize--
	}

	// create the FIR Parks-McClellan struct
	firpm = newFirpm(gridsize, r)

	/*
	* Create dense frequency grid
	 */
	firpm.createDenseGrid(r, numtaps, gridsize, density, f, a, w, symmetry)
	firpm.initialGuess(r, gridsize)

	/*
	* For Differentiator: (fix grid)
	 */
	if ftype == differentiator {
		for i := 0; i < gridsize; i++ {
			/* D[i] = D[i]*Grid[i]; */
			if firpm.D[i] > 0.0001 {
				firpm.W[i] = firpm.W[i] / firpm.Grid[i]
			}
		}
	}

	/*
	* For odd or negative symmetry filters, alter the
	* D[] and W[] according to Parks McClellan
	 */
	var c float64
	if symmetry == positive {
		if numtaps%2 == 0 {
			for i := 0; i < gridsize; i++ {
				c = math.Cos(pi * firpm.Grid[i])
				firpm.D[i] /= c
				firpm.W[i] *= c
			}
		}
	} else {
		if numtaps%2 != 0 {
			for i := 0; i < gridsize; i++ {
				c = math.Sin(twoPi * firpm.Grid[i])
				firpm.D[i] /= c
				firpm.W[i] *= c
			}
		} else {
			for i := 0; i < gridsize; i++ {
				c = math.Sin(pi * firpm.Grid[i])
				firpm.D[i] /= c
				firpm.W[i] *= c
			}
		}
	}

	/*
	* Perform the Remez Exchange algorithm
	 */
	for iter := 0; iter < maxIterations; iter++ {
		firpm.calcParams(r)
		firpm.calcError(r, gridsize)
		err := firpm.search(r, gridsize)
		if err != nil {
			return nil, err
		}

		if firpm.isDone(r) {
			break
		}

	}

	firpm.calcParams(r)

	/*
	* Find the 'taps' of the filter for use with Frequency
	* Sampling.  If odd or Negative symmetry, fix the taps
	* according to Parks McClellan
	 */
	for i := 0; i <= numtaps/2; i++ {
		if symmetry == positive {
			if numtaps%2 != 0 {
				c = 1.0
			} else {
				c = math.Cos(pi * float64(i) / float64(numtaps))
			}
		} else {
			if numtaps%2 != 0 {
				c = math.Sin(twoPi * float64(i) / float64(numtaps))
			} else {
				c = math.Sin(pi * float64(i) / float64(numtaps))
			}
		}
		firpm.taps[i] = firpm.computeA(float64(i)/float64(numtaps), r) * c
	}

	/*
	* Frequency sampling design with calculated taps
	 */
	firpm.freqSample(numtaps, firpm.taps, h, symmetry)

	return h, nil
}
