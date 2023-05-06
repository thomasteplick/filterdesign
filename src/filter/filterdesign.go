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
	"bufio"
	"fmt"
	"html/template"
	"log"
	"math"
	"math/cmplx"
	"net/http"
	"os"
	"path"
	"strconv"

	"github.com/mjibson/go-dsp/fft"
	"github.com/thomasteplick/firpm"
)

const (
	rows                       = 300                                  // #rows in grid
	columns                    = 300                                  // #columns in grid
	addr                       = "127.0.0.1:8080"                     // http server listen address
	filterdesignoptions        = "templates/filterdesignoptions.html" // html for filter design options
	plotfilter                 = "templates/plotfilter.html"          // html for filter plotting in time or frequency domain
	patternFilterDesignOptions = "/filterdesignoptions"               // http handler for filter design options
	patternPlotFilter          = "/plotfilter"                        // http handler for filter plot in time or frequency domain
	dataDir                    = "data/"                              // directory for the data files
	xlabels                    = 11                                   // # labels on x axis
	ylabels                    = 11                                   // # labels on y axis
	deg2rad                    = math.Pi / 180.0                      // convert degrees to radians
	FFTSize                    = 8192                                 // default FFT size
)

// filter types
const (
	lowpass        = iota // 0, not used
	bandpass              // 1
	differentiator        // 2
	hilbert               // 3
	highpass              // 4, not used
)

// band properties
const (
	freq1  = iota // 0 start edge
	freq2         // 1 end edge
	ampl1         // 2 amplitude at start edge
	ampl2         // 3 amplitude at end edge
	weight        // 4 weight of the edge
)

// Type to contain all the HTML template actions
type PlotT struct {
	Grid         []string // plotting grid
	Status       string   // status of the plot
	Xlabel       []string // x-axis labels
	Ylabel       []string // y-axis labels
	Filename     string   // filename to plot
	SamplingFreq string   // data sampling rate in Hz
	FFTSegments  string   // FFT segments, K
	FFTSize      string   // FFT size
	Samples      string   // complex samples in data file
}

// Type to hold the minimum and maximum data values
type Endpoints struct {
	xmin float64
	xmax float64
	ymin float64
	ymax float64
}

// band property attributes
type Attributes struct {
	Name  string
	Value string
	Class string
}

// frequency band properties
type Band []Attributes

// HTML form controls to be executed with HTML template
type FormControls struct {
	Status      string // form messages resulting from a submit
	State       string // valid or invalid, depending on band properties
	FilterOrder string // > 2
	filterType  string // lowpass, highpass, bandpass, hilbert, differentiator
	numbands    string // > 1
	Density     string // default 16
	Bands       []Band // <= 5 filter bands, table of band properties
}

// global variables for parse and execution of the html template
var (
	tmplFormFilterDesignOptions *template.Template
	tmplPlotFilter              *template.Template
	bandProperties              = []int{freq1, freq2, ampl1, ampl2, weight}
)

// init parses the HTML template file
func init() {
	tmplFormFilterDesignOptions = template.Must(template.ParseFiles(filterdesignoptions))
	tmplPlotFilter = template.Must(template.ParseFiles(plotfilter))
}

// handleFilterDesignOptions creates the FIR filter using Parks-McClellan algorithm
func handleFilterDesignOptions(w http.ResponseWriter, r *http.Request) {
	/*
			Tips
			If your filter design fails to converge, the filter design might not be correct.
			Verify the design by checking the frequency response.

			If your filter design fails to converge and the resulting filter design is not correct,
		     attempt one or more of the following:

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
	if len(temp) > 0 {
		// verify filter order is > 2
		order, err := strconv.Atoi(temp)
		form.FilterOrder = temp
		if err != nil {
			form.Status = fmt.Sprintf("Filter order conversion error: %v", err)
		} else {
			if order < 3 {
				form.Status = "Filter order must be greater than 2"
			}
		}

		// Get filter type, change lowpass or highpass to bandpass
		// verify ftype in ["bandpass", "hilbert", "differentiator"]
		filterType := r.FormValue("filtertype")
		if len(filterType) == 0 {
			form.Status = "Filter type must be specified"
		} else {
			if filterType == "highpass" || filterType == "lowpass" ||
				filterType == "bandpass" || filterType == "multibandpass" {
				ftype = bandpass
			} else if filterType == "differentiator" {
				ftype = differentiator
			} else if filterType == "hilbert" {
				ftype = hilbert
			} else {
				form.Status = fmt.Sprintf("filter type %s not allowed", filterType)
			}
		}
		form.filterType = filterType

		// Get density and verify >= 16
		temp = r.FormValue("griddensity")
		if len(temp) == 0 {
			form.Status = "grid density not defined"
		} else {
			density, err = strconv.Atoi(temp)
			if err != nil {
				form.Status = fmt.Sprintf("grid density int conversion error: %v", err)
			} else {
				form.Density = temp
				if density < 16 {
					form.Status = "grid density must be >= 16"
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
			form.Status = fmt.Sprintf("Number of bands conversion error: %v", err)
		} else if numband < 2 {
			form.Status = "number of frequency bands must be greater than one"
		} else {
			form.Bands = make([]Band, numband)
			for i := range form.Bands {
				form.Bands[i] = make([]Attributes, len(bandProperties))
			}
			// slices for Remez call parameters
			freq := make([]float64, 2*numband)
			ampl := make([]float64, 2*numband)
			wt := make([]float64, numband)
			// loop over the bands, get and verify properties, and store in the form struct
			for n := 0; n < numband; n++ {
				// loop over the band properties, get band weights from table
				// verify len(w) == len(f)/2
				for _, prop := range bandProperties {
					name := fmt.Sprintf("band%d_prop%d", n, prop)
					temp = r.FormValue(name)
					if len(temp) == 0 {
						form.Status = fmt.Sprintf("missing band %d property %d", n, prop)
						form.Bands[n][prop] = Attributes{Name: name, Class: "invalid"}
					} else {
						param, err := strconv.ParseFloat(temp, 64)
						if err != nil {
							form.Status = fmt.Sprintf("band %d property %d float conversion error: %v", n, prop, err)
							form.Bands[n][prop] = Attributes{Name: name, Class: "invalid"}
						} else {
							if prop == freq1 {
								form.Bands[n][freq1] = Attributes{Name: name, Class: "valid", Value: temp}
								// firpm.Remez requires 0 <= normalized frequency <= 0.5
								freq[2*n] = param / 2.0
							} else if prop == freq2 {
								form.Bands[n][freq2] = Attributes{Name: name, Class: "valid", Value: temp}
								// firpm.Remez requires 0 <= normalized frequency <= 0.5
								freq[2*n+1] = param / 2.0
							} else if prop == ampl1 {
								form.Bands[n][ampl1] = Attributes{Name: name, Class: "valid", Value: temp}
								ampl[2*n] = param
							} else if prop == ampl2 {
								form.Bands[n][ampl2] = Attributes{Name: name, Class: "valid", Value: temp}
								ampl[2*n+1] = param
							} else {
								form.Bands[n][weight] = Attributes{Name: name, Class: "valid", Value: temp}
								wt[n] = param
							}
						}
					}
				}
			}
			// verify 0 <= freq[i] <= 1, normalized frequency where 1 = Nyquist frequency
			// verify the first frequency is 0 and the last 1.0
			for n := 0; n < 2*numband-1; n += 2 {
				if freq[n] < 0 || freq[n] > 0.5 {
					form.Status = fmt.Sprintf("band %d start edge frequency not >=0 and <= 1.0", n/2)
					form.Bands[n/2][freq1].Class = "invalid"
				}
				if freq[n+1] < 0 || freq[n+1] > 0.5 {
					form.Status = fmt.Sprintf("band %d end edge frequency not >=0 and <= 1.0", n/2)
					form.Bands[n/2][freq2].Class = "invalid"
				}
			}
			if freq[0] != 0.0 {
				form.Status = "First start edge frequency must be 0"
				form.Bands[0][freq1].Class = "invalid"
			}
			if freq[2*numband-1] != 0.5 {
				form.Status = "Last end edge frequency must be 1.0"
				form.Bands[numband-1][freq2].Class = "invalid"
			}

			// verify freq[i] < freq[i+1]
			for n := 0; n < 2*numband-1; n++ {
				if freq[n] > freq[n+1] {
					if n%2 == 0 {
						form.Bands[n/2][freq1].Class = "invalid"
						form.Bands[n/2][freq2].Class = "invalid"
						form.Status = fmt.Sprintf("band %d freq1 and freq2 not in the right order", n/2)
					} else {
						form.Bands[n/2][freq2].Class = "invalid"
						form.Bands[(n+1)/2][freq1].Class = "invalid"
						form.Status = fmt.Sprintf("band %d freq2 and band %d freq1 not in the right order", n/2, (n+1)/2)
					}
				}
			}

			// Get Remez filter coefficients
			if len(form.Status) == 0 {
				h, err = firpm.Remez(order, freq, ampl, wt, ftype, density)
				if err != nil {
					// update status with error
					form.Status = err.Error()
				} else {
					// update status with coefficient file location
					// Save the FIR coefficients
					f, err := os.Create(path.Join(dataDir, filename))
					if err != nil {
						form.Status = fmt.Sprintf("create %v error: %v", path.Join(dataDir, filename), err)
					} else {
						defer f.Close()
						for i := 0; i <= order; i++ {
							fmt.Fprintf(f, "%v\n", h[i])
						}
					}
				}
			}
			if len(form.Status) > 0 {
				form.State = "invalid"
			} else {
				form.State = "valid"
				form.Status = fmt.Sprintf("FIR coefficients saved in %s", filename)
			}
			// Write to HTTP using template and form
			if err := tmplFormFilterDesignOptions.Execute(w, form); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
		}
	} else {
		form.Status = "Fill in the required filter parameters and submit"
		// Write to HTTP using template and form``
		if err := tmplFormFilterDesignOptions.Execute(w, form); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}
	}
}

// findEndpoints finds the minimum and maximum filter values
func (ep *Endpoints) findEndpoints(input *bufio.Scanner) {
	ep.xmax = 0.0
	ep.xmin = 0.0
	ep.ymax = -math.MaxFloat64
	ep.ymin = math.MaxFloat64
	for input.Scan() {
		value := input.Text()
		var (
			y   float64
			err error
		)

		if y, err = strconv.ParseFloat(value, 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", value, err)
			continue
		}

		if y > ep.ymax {
			ep.ymax = y
		}
		if y < ep.ymin {
			ep.ymin = y
		}

		ep.xmax++
	}
	ep.xmax--
}

// gridFillInterp inserts the data points in the grid and draws a straight line between points
func gridFillInterp(plot *PlotT, xscale float64, yscale float64, endpoints Endpoints, input *bufio.Scanner) error {

	var (
		x            float64 = 0.0
		y            float64 = 0.0
		prevX, prevY float64
		err          error
	)

	// Get first sample
	input.Scan()
	value := input.Text()

	if y, err = strconv.ParseFloat(value, 64); err != nil {
		fmt.Printf("String %s conversion to float error: %v\n", value, err)
		return err
	}

	// This cell location (row,col) is on the line
	row := int((endpoints.ymax-y)*yscale + .5)
	col := int((x-endpoints.xmin)*xscale + .5)
	plot.Grid[row*columns+col] = "online"

	prevX = x
	prevY = y

	// Scale factor to determine the number of interpolation points
	/*
		lenEP := math.Sqrt((endpoints.xmax-endpoints.xmin)*(endpoints.xmax-endpoints.xmin) +
			(endpoints.ymax-endpoints.ymin)*(endpoints.ymax-endpoints.ymin))
	*/
	lenEPy := endpoints.ymax - endpoints.ymin
	lenEPx := endpoints.xmax - endpoints.xmin

	// Continue with the rest of the points in the file
	for input.Scan() {
		x++
		value = input.Text()
		if y, err = strconv.ParseFloat(value, 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", value, err)
			return err
		}

		// This cell location (row,col) is on the line
		row := int((endpoints.ymax-y)*yscale + .5)
		col := int((x-endpoints.xmin)*xscale + .5)
		plot.Grid[row*columns+col] = "online"

		// Interpolate the points between previous point and current point

		/* lenEdge := math.Sqrt((x-prevX)*(x-prevX) + (y-prevY)*(y-prevY)) */
		lenEdgeX := math.Abs((x - prevX))
		lenEdgeY := math.Abs(y - prevY)
		ncellsX := int(columns * lenEdgeX / lenEPx) // number of points to interpolate in x-dim
		ncellsY := int(rows * lenEdgeY / lenEPy)    // number of points to interpolate in y-dim
		// Choose the biggest
		ncells := ncellsX
		if ncellsY > ncells {
			ncells = ncellsY
		}

		stepX := (x - prevX) / float64(ncells)
		stepY := (y - prevY) / float64(ncells)

		// loop to draw the points
		interpX := prevX
		interpY := prevY
		for i := 0; i < ncells; i++ {
			row := int((endpoints.ymax-interpY)*yscale + .5)
			col := int((interpX-endpoints.xmin)*xscale + .5)
			plot.Grid[row*columns+col] = "online"
			interpX += stepX
			interpY += stepY
		}

		// Update the previous point with the current point
		prevX = x
		prevY = y
	}
	return nil
}

// processTimeDomain plots the time domain data from disk file
func processTimeDomain(w http.ResponseWriter, r *http.Request, filename string) error {

	// main data structure
	var (
		plot      PlotT
		xscale    float64
		yscale    float64
		endpoints Endpoints
	)

	plot.Grid = make([]string, rows*columns)
	plot.Xlabel = make([]string, xlabels)
	plot.Ylabel = make([]string, ylabels)

	// Open file
	f, err := os.Open(filename)
	if err == nil {
		// Mark the data x-y coordinate online at the corresponding
		// grid row/column.
		input := bufio.NewScanner(f)

		endpoints.findEndpoints(input)

		f.Close()
		f, err = os.Open(filename)
		if err == nil {
			defer f.Close()
			input := bufio.NewScanner(f)

			// Calculate scale factors for x and y
			xscale = (columns - 1) / (endpoints.xmax - endpoints.xmin)
			yscale = (rows - 1) / (endpoints.ymax - endpoints.ymin)

			// Fill in the grid with the data points using interpolation
			err = gridFillInterp(&plot, xscale, yscale, endpoints, input)
			if err != nil {
				return err
			}

			// Set plot status if no errors
			if len(plot.Status) == 0 {
				plot.Status = fmt.Sprintf("Status: Data plotted from (%.3f,%.3f) to (%.3f,%.3f)",
					endpoints.xmin, endpoints.ymin, endpoints.xmax, endpoints.ymax)
			}

		} else {
			// Set plot status
			fmt.Printf("Error opening file %s: %v\n", filename, err)
			return fmt.Errorf("error opening file %s: %v", filename, err)
		}
	} else {
		// Set plot status
		fmt.Printf("Error opening file %s: %v\n", filename, err)
		return fmt.Errorf("error opening file %s: %v", filename, err)
	}

	// Construct x-axis labels
	incr := (endpoints.xmax - endpoints.xmin) / (xlabels - 1)
	x := endpoints.xmin
	// First label is empty for alignment purposes
	for i := range plot.Xlabel {
		plot.Xlabel[i] = fmt.Sprintf("%.2f", x)
		x += incr
	}

	// Construct the y-axis labels
	incr = (endpoints.ymax - endpoints.ymin) / (ylabels - 1)
	y := endpoints.ymin
	for i := range plot.Ylabel {
		plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
		y += incr
	}

	// Enter the filename in the form
	plot.Filename = path.Base(filename)

	// Write to HTTP using template and grid
	if err := tmplPlotFilter.Execute(w, plot); err != nil {
		log.Fatalf("Write to HTTP output using template with error: %v\n", err)
	}

	return nil
}

// processFrequencyDomain calculates the power spectral density (PSD) and plots it
func processFrequencyDomain(w http.ResponseWriter, r *http.Request, filename string) error {
	// Use complex128 for FFT computation

	var (
		plot         PlotT // main data structure to execute with parsed html template
		endpoints    Endpoints
		N            int                          //  complex FFT size
		nn           int                          // number of complex samples in the data file
		K            int                          //  number of segments used in PSD with 50% overlap
		PSD          []float64                    // power spectral density
		psdMax       float64   = -math.MaxFloat64 // maximum PSD value
		psdMin       float64   = math.MaxFloat64  // minimum PSD value
		xscale       float64                      // data to grid in x direction
		yscale       float64                      // data to grid in y direction
		samplingRate float64   = 2.0              // sampling rate in Hz
	)

	plot.Grid = make([]string, rows*columns)
	plot.Xlabel = make([]string, xlabels)
	plot.Ylabel = make([]string, ylabels)

	N = FFTSize

	// Power Spectral Density, PSD[N/2] is the Nyquist critical frequency
	// It is (sampling frequency)/2, the highest non-aliased frequency
	PSD = make([]float64, N/2)

	// Open the filter file
	f, err := os.Open(filename)
	if err == nil {
		defer f.Close()
		bufN := make([]complex128, N)
		input := bufio.NewScanner(f)
		// Read in all the samples
		k := 0
		var real float64
		for input.Scan() {
			value := input.Text()
			if real, err = strconv.ParseFloat(value, 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", value, err)
				continue
			}

			bufN[k] = complex(real, 0)
			k++
		}

		// zero-pad N-k samples in bufN
		for i := k; i < N; i++ {
			bufN[i] = 0
		}

		// Perform N-point complex FFT and insert magnitude^2 in PSD
		// Then convert to dB with 10*log10()
		fourierN := fft.FFT(bufN)
		x := cmplx.Abs(fourierN[0])
		PSD[0] = 10.0 * math.Log10(x*x)
		psdMax = PSD[0]
		psdMin = PSD[0]
		for i := 1; i < N/2; i++ {
			// Use positive and negative frequencies -> bufN[N-i] = bufN[-i]
			xi := cmplx.Abs(fourierN[i])
			xNi := cmplx.Abs(fourierN[N-i])
			PSD[i] = 10.0 * math.Log10(xi*xi+xNi*xNi)
			if PSD[i] > psdMax {
				psdMax = PSD[i]
			}
			if PSD[i] < psdMin {
				psdMin = PSD[i]
			}
		}

		endpoints.xmin = 0
		endpoints.xmax = float64(N / 2) // equivalent to Nyquist critical frequency
		endpoints.ymin = psdMin
		endpoints.ymax = psdMax

		// Calculate scale factors for x and y
		xscale = (columns - 1) / (endpoints.xmax - endpoints.xmin)
		yscale = (rows - 1) / (endpoints.ymax - endpoints.ymin)

		// Store the PSD in the plot Grid
		for bin, pow := range PSD {
			row := int((endpoints.ymax-pow)*yscale + .5)
			col := int((float64(bin)-endpoints.xmin)*xscale + .5)
			plot.Grid[row*columns+col] = "online"
		}

		// Plot the PSD N/2 float64 values, execute the data on the plotfrequency.html template

		// Set plot status if no errors
		if len(plot.Status) == 0 {
			plot.Status = fmt.Sprintf("Status: Data plotted from (%.3f,%.3f) to (%.3f,%.3f)",
				endpoints.xmin, endpoints.ymin, endpoints.xmax, endpoints.ymax)
		}

	} else {
		// Set plot status
		fmt.Printf("Error opening file %s: %v\n", filename, err)
		return fmt.Errorf("error opening file %s: %v", filename, err)
	}

	// Apply the  sampling rate in Hz to the x-axis using a scale factor
	// Convert the fft size to sampling rate/2, the Nyquist critical frequency
	sf := 0.5 * samplingRate / endpoints.xmax

	// Construct x-axis labels
	incr := (endpoints.xmax - endpoints.xmin) / (xlabels - 1)
	x := endpoints.xmin
	// First label is empty for alignment purposes
	for i := range plot.Xlabel {
		plot.Xlabel[i] = fmt.Sprintf("%.2f", x*sf)
		x += incr
	}

	// Construct the y-axis labels
	incr = (endpoints.ymax - endpoints.ymin) / (ylabels - 1)
	y := endpoints.ymin
	for i := range plot.Ylabel {
		plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
		y += incr
	}

	// Insert frequency domain parameters in the form
	plot.SamplingFreq = fmt.Sprintf("%.0f", samplingRate)
	plot.FFTSegments = strconv.Itoa(K)
	plot.FFTSize = strconv.Itoa(N)
	plot.Samples = strconv.Itoa(nn)

	// Enter the filename in the form
	plot.Filename = path.Base(filename)

	// Write to HTTP using template and grid
	if err := tmplPlotFilter.Execute(w, plot); err != nil {
		log.Fatalf("Write to HTTP output using template with error: %v\n", err)
	}

	return nil
}

// handlePlotFilter plots the filter in time or frequency domains
func handlePlotFilter(w http.ResponseWriter, r *http.Request) {
	// main data structure
	var (
		plot PlotT
		err  error = nil
	)

	filename := r.FormValue("filename")
	// choose time or frequency domain processing
	if len(filename) > 0 {

		domain := r.FormValue("domain")
		switch domain {
		case "time":
			err = processTimeDomain(w, r, path.Join(dataDir, filename))
			if err != nil {
				plot.Status = err.Error()
			}
		case "frequency":
			err = processFrequencyDomain(w, r, path.Join(dataDir, filename))
			if err != nil {
				plot.Status = err.Error()
			}
		default:
			plot.Status = fmt.Sprintf("Invalid domain choice: %s", domain)
		}

		if err != nil {

			// Write to HTTP using template and grid
			if err := tmplPlotFilter.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
		}
	} else {
		plot.Status = "Enter filename and domain to plot"
		// Write to HTTP using template and grid``
		if err := tmplPlotFilter.Execute(w, plot); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}
	}
}

// main sets up the http handlers, listens, and serves http clients
func main() {
	// Set up http servers with handler for Filter Design Options and Plot Filter
	http.HandleFunc(patternFilterDesignOptions, handleFilterDesignOptions)
	http.HandleFunc(patternPlotFilter, handlePlotFilter)
	fmt.Printf("Filter Design Server listening on %v.\n", addr)
	http.ListenAndServe(addr, nil)
}
