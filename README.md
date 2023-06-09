# filterdesign
Design a FIR filter using the Parks-McClellan algorithm

This is a web application written in Go which uses the html/template package.  It can be accessed at http://127.0.0.1:8080/filterdesignoptions.
The user enters the filter order, grid density, filter type, and number of bands.  After clicking on the
submit button, a frequency band table appears. Fill in the begin and end frequency of each band, the amplitudes
at the edges of each band, and the weights given for each band.  The frequencies are normalized
between 0 and 1, with 1 being the Nyquist frequency (1/2 the sampling frequency).  The weight
gives the ability to specify the error reduction between the desired and actual amplitudes
for each band.  The bigger the number, the smaller the error between the desired and actual 
amplitudes.  The frequencies must start at 0 and end at 1, in increasing order.  The impulse and frequency
responses can be displayed by clicking on the Plot Filter link.  Choose the file and time or frequency domain and enter submit.
The frequency response is in dB.  The user can create another design by clicking on the 
Filter Design Options link.  The Grid Density determines how accurately the filter will be constructed.
The minimum value is 16, but higher numbers are slower to compute.  The maximum number of bands that can be designed 
is 5, and the minimum is 2.  The filter coefficients are stored in the data folder.  The filter types are
<ul>
<li>Lowpass</li>
<li>Highpass</li>
<li>Bandpass</li>
<li>Multibandpass</li>
<li>Differentiator</li>
<li>Hilbert Transformer</li>
</ul>
The Parks-McClellan algorithm uses the Remez exchange algorithm and Chebyshev
approximation theory to design filters with optimal fits between the
desired and actual frequency responses. The filters are optimal in
the sense that the maximum error between the desired frequency response
and the actual frequency response is minimized.

<h3> Filter Desigin Options</h3>

![image](https://user-images.githubusercontent.com/117768679/236644927-ff193723-b18d-4437-9d30-905a4af8e204.png)

<h3> Frequency Band Table Before Submission</h3>

![image](https://user-images.githubusercontent.com/117768679/236645232-2cde0e19-8877-4eb4-8614-7fe740440336.png)

<h3>Frequency Band Table After Submission</h3>

![image](https://user-images.githubusercontent.com/117768679/236645427-c0f39dec-7e2f-4ad6-b015-f6b62f8ec4ab.png)

<h3>Impulse Response</h3>>
  
![image](https://user-images.githubusercontent.com/117768679/236645585-a0e089d6-3a2d-414d-ab0a-efff21404f68.png)
  
<h3>Frequency Response</h3>

![image](https://user-images.githubusercontent.com/117768679/236645784-37288bb4-e66f-40e5-87dc-19d39ceb00c1.png)
  
 <h3>Multibandpass Filter Before Submit</h3>
  
 ![image](https://user-images.githubusercontent.com/117768679/236700118-eef1d80a-5aaa-407f-af5b-605fea80ae89.png)
  
 <h3>Multibandpass Filter After Submit</h3>
  
 ![image](https://user-images.githubusercontent.com/117768679/236700219-bcb12bc8-24f6-43d1-99b6-382c1f091ad9.png)

 <h3>Impulse Response</h3>
  
![image](https://user-images.githubusercontent.com/117768679/236700469-b6f87c71-b7c2-4e25-bd5e-0ba1bf394ad7.png)
  
<h3>Frequency Response</>
  
![image](https://user-images.githubusercontent.com/117768679/236700587-398563e8-ed96-4594-87cc-180392576b4f.png)

<h3>Low Pass Filter before submission</h3>

![image](https://user-images.githubusercontent.com/117768679/236735377-4187269f-01f3-41ce-ac3f-ac29d1e9ee3e.png)

<h3>Low Pass Filter after submission</h3>

![image](https://user-images.githubusercontent.com/117768679/236735723-7aeb48e2-4e4a-4ac0-b70e-246774e84e88.png)

<h3>Impulse Response</h3>

![image](https://user-images.githubusercontent.com/117768679/236735061-5f82d551-e79a-4b6a-9faa-e3f739b38b6a.png)

<h3>Filter Response</h3>

![image](https://user-images.githubusercontent.com/117768679/236734905-83ee92d3-9777-404e-ad29-83a2ebc955f2.png)

<h3>Multibandpass Filter Before Submit</h3>

![image](https://github.com/thomasteplick/filterdesign/assets/117768679/c9fee99f-2d89-48e9-85d5-a747faadb90b)

<h3>Multibandpass Filter After Submit</h3>

![image](https://github.com/thomasteplick/filterdesign/assets/117768679/91a6d205-b8ef-4f8d-8b54-e9d234bca035)

<h3>Impulse Response</>
  
![image](https://github.com/thomasteplick/filterdesign/assets/117768679/23866286-4a78-4c55-b015-22295ce938ae)

<h3>Filter Response</h3>

![image](https://github.com/thomasteplick/filterdesign/assets/117768679/41720a72-da0a-4479-8572-c518cfa7ad00)

