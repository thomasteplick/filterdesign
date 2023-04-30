/*
 * filterdesign
 * Create a filter according to user's choices.  A Web application in which
 * the user enters the filter specs and submits them in an HTML form.
 * Currently supported are:
 *    Parks-McClellan algorithm for equiripple FIR filter.  The user enters
 *    the filter order, type, normalized frequences, amplitudes in the bands,
 *    and band weights.
 *
 *  Parks-McClellan algorithm uses the Remez exchange algorithm and Chebyshev
 *  approximation theory to design filters with optimal fits between the
 *  desired and actual frequency responses. The filters are optimal in
 *  the sense that the maximum error between the desired frequency response
 *  and the actual frequency response is minimized.
 */

package main

import (
	"fmt"
	"html/template"
	"log"
	"net/http"
	"os"
	"path"
	"strconv"

	"github.com/thomasteplick/firpm"
)

const (
	addr                       = "127.0.0.1:8080"                     // http server listen address
	filterdesignoptions        = "templates/filterdesignoptions.html" // html for filter design options
	patternFilterDesignOptions = "/filterdesignoptions"               // http handler for filter design options
	dataDir                    = "data/"                              // directory for the data files
)

// filter types
const (
	lowpass        = iota // 0, not used
	bandpass              // 1
	differentiator        // 2
	hilbert               // 3
	highpass              // 4, not used
)

// frequency band properties
type Band struct {
	state string // valid or invalid
	freq1 string // start edge frequency (normalized) of band
	freq2 string // end edge frequency (normalized) of band
	ampl1 string // start edge amplitude
	ampl2 string // end edge amplitude
	wt    string // band weight
}

// HTML form controls to be executed with HTML template
type FormControls struct {
	status      string // form messages resulting from a submit
	state       string // valid or invalid, depending on band properties
	filterOrder string // > 2
	filterType  string // lowpass, highpass, bandpass, hilbert, differentiator
	numbands    string // > 1
	density     string // default 16
	bands       []Band // <= 5 filter bands, table of band properties
}

// global variables for parse and execution of the html template
var (
	tmplFormFilterDesignOptions *template.Template
	bandProperties              = []string{"f1", "f2", "a1", "a2", "wt"}
)

// init parses the HTML template file
func init() {
	tmplFormFilterDesignOptions = template.Must(template.ParseFiles(filterdesignoptions))
}

// handleFilterDesignOptions
func handleFilterDesignOptions(w http.ResponseWriter, r *http.Request) {
	/*
		Tips
		If your filter design fails to converge, the filter design might not be correct. Verify the design by checking the frequency response.

		If your filter design fails to converge and the resulting filter design is not correct, attempt one or more of the following:

		Increase the filter order.

		Relax the filter design by reducing the attenuation in the stopbands and/or broadening the transition regions.
	*/

	var (
		form     FormControls // to be executed on HTML template
		h        []float64    // FIR filter coefficients
		filename string       // filter coefficients storage
		density  int          // grid density >= 16
		ftype    int          // filter type integer
	)

	// Get filter order
	temp := r.FormValue("filterorder")
	form.filterOrder = temp
	if len(temp) > 0 {
		// verify filter order is > 2
		order, err := strconv.Atoi(temp)
		if err != nil {
			form.status = fmt.Sprintf("Filter order conversion error: %v", err)
		} else {
			if order < 3 {
				form.status = "Filter order must be greater than 2"
			}
		}

		// Get filter type, change lowpass or highpass to bandpass
		// verify ftype in ["bandpass", "hilbert", "differentiator"]
		filterType := r.FormValue("filtertype")
		form.filterType = filterType
		if len(filterType) == 0 {
			form.status = "Filter type must be specified"
		} else {
			if filterType == "highpass" || filterType == "lowpass" ||
				filterType == "bandpass" {
				filterType = "bandpass"
				ftype = bandpass
			} else if filterType == "differentiator" {
				ftype = differentiator
			} else if filterType == "hilbert" {
				ftype = hilbert
			} else {
				form.status = fmt.Sprintf("filter type %s not allowed", filterType)
			}
		}

		// Get density and verify >= 16
		temp = r.FormValue("griddensity")
		if len(temp) == 0 {
			form.status = "grid density not defined"
		} else {
			density, err = strconv.Atoi(temp)
			if err != nil {
				form.status = fmt.Sprintf("grid density int conversion error: %v", err)
			} else {
				form.density = temp
				if density < 16 {
					form.status = "grid density must be >= 16"
				}
			}
		}

		// firpm{filterType}{order}.txt
		filename = fmt.Sprintf("firpm%s%d.txt", filterType, order)

		// Get frequency bands from table and verify numbands > 1
		temp = r.FormValue("numbands")
		form.numbands = temp
		numband, err := strconv.Atoi(temp)
		if err != nil {
			form.status = fmt.Sprintf("Number of bands conversion error: %v", err)
		} else if numband < 2 {
			form.status = "number of frequency bands must be greater than one"
		} else {
			form.bands = make([]Band, numband)
			// slices for Remez call parameters
			freq := make([]float64, 2*numband)
			ampl := make([]float64, 2*numband)
			wt := make([]float64, numband)
			// loop over the bands, get and verify properties, and store in the form struct
			for n := 0; n < numband; n++ {
				// loop over the band properties, get band weights from table
				// verify len(w) == len(f)/2
				for _, prop := range bandProperties {
					temp = r.FormValue("band" + strconv.Itoa(n) + prop)
					if len(temp) == 0 {
						form.status = fmt.Sprintf("missing band %d property %s", n+1, prop)
						form.bands[n].state = "invalid"
					} else {
						param, err := strconv.ParseFloat(temp, 64)
						if err != nil {
							form.status = fmt.Sprintf("band %d property %s float conversion error: %v", n+1, prop, err)
							form.bands[n].state = "invalid"
						} else {
							if prop == "f1" {
								form.bands[n].freq1 = temp
								freq[2*n] = param
							} else if prop == "f2" {
								form.bands[n].freq2 = temp
								freq[2*n+1] = param
							} else if prop == "a1" {
								form.bands[n].ampl1 = temp
								ampl[2*n] = param
							} else if prop == "a2" {
								form.bands[n].ampl2 = temp
								ampl[2*n+1] = param
							} else {
								form.bands[n].wt = temp
								wt[n] = param
							}
						}
					}
				}
			}
			// verify 0 <= freq[i] <= 1, normalized frequency where 1 = Nyquist frequency
			// verify the first frequency is 0 and the last 1
			for n := 0; n < 2*numband; n += 2 {
				if freq[n] < 0 || freq[n] > 1.0 {
					form.status = fmt.Sprintf("band %d start edge frequency not >=0 and <= 1.0", n/2+1)
					form.bands[n].state = "invalid"
				}
				if freq[n+1] < 0 || freq[n+1] > 1.0 {
					form.status = fmt.Sprintf("band %d end edge frequency not >=0 and <= 1.0", n/2+1)
					form.bands[n+1].state = "invalid"
				}
			}
			if freq[0] != 0.0 {
				form.status = "First start edge frequency must be 0"
				form.bands[0].state = "invalid"
			}
			if freq[2*numband-1] != 1.0 {
				form.status = "Last end edge frequency must be 1.0"
				form.bands[numband-1].state = "invalid"
			}

			// verify freq[i] < freq[i+1]
			for i := 0; i < 2*numband-1; i++ {
				if freq[i] > freq[i+1] {
					form.bands[i/2].state = "invalid"
					form.status = fmt.Sprintf("band %d frequency > band frequency %d", i/2, (i+1)/2)
				}
			}

			// Get Remez filter coefficients
			if len(form.status) == 0 {
				h, err = firpm.Remez(order, freq, ampl, wt, ftype, density)
				if err != nil {
					// update status with error
					form.status = err.Error()
				} else {
					// update status with coefficient file location
					form.status = fmt.Sprintf("FIR coefficients saved in %s", filename)
					// Save the FIR coefficients
					f, err := os.Create(path.Join(dataDir, filename))
					if err != nil {
						form.status = fmt.Sprintf("create %v error: %v", path.Join(dataDir, filename), err)
					} else {
						defer f.Close()
						// write the time so that the plotting program can figure out the sampling frequency
						// this will allow it to label the x-axis with the normalized frequency 0 -> 1
						t := 0.0
						samplingRate := 2.0 // in Hz, gives a Nyquist frequency of 1 Hz
						step := 1.0 / samplingRate
						for i := 0; i <= order; i++ {
							fmt.Fprintf(f, "%v %v\n", t, h[i])
							t += step
						}
					}
				}
			}
			form.state = "valid"
			if len(form.status) > 0 {
				form.state = "invalid"
			}
			// Write to HTTP using template and form
			if err := tmplFormFilterDesignOptions.Execute(w, form); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
		}
	} else {
		form.status = "Fill in the required filter parameters and submit"
		// Write to HTTP using template and form``
		if err := tmplFormFilterDesignOptions.Execute(w, form); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}
	}
}

// main sets up the http handlers, listens, and serves http clients
func main() {
	// Set up http servers with handler for Filter Design Options
	http.HandleFunc(patternFilterDesignOptions, handleFilterDesignOptions)
	fmt.Printf("Filter Design Server listening on %v.\n", addr)
	http.ListenAndServe(addr, nil)
}
