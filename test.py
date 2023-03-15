#! /usr/bin/python
import ctypes
import glob
import numpy
import scipy.fftpack #Alternative FFT
import matplotlib.pyplot as plt
import sys
import time


from operator import add
from matplotlib import pyplot
from matplotlib import mlab
from ctypes import *
from numpy.ctypeslib import ndpointer
from math import log
from sklearn.linear_model import LinearRegression

steps = 2**17 #A singal composed of n steps
fhigh = 1000  #Sampling frequency, Nyquist should be = 0.5*fhigh
flow = 0.1    #Will determine number of decades, e.g. log10(fhigh/flow) 
Gain = 1      #Scalar for the noise
h = 5         #Numbr of poles, >>h = >>filter coeffecients 

libfile = glob.glob('build/*/noise*.so')[0]

mylib = ctypes.CDLL(libfile)

mylib.noise.restype = ndpointer(dtype=ctypes.c_double, shape=(steps,))
mylib.noise.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]

num_to_avg = 1000

start_time = time.time()

#This loop goes through j=3 for gamma=0,1,2 or j=n/b for real values of gamma. 
#Then produces num_to_avg signals based on each gamma.
#Each averaged signal is later plotted and/or its psd

for j in range(5):
	
	steps_2 = int((steps/2))
	psd_sum = [0]*steps_2
	noise_sum = [0]*steps

	for i in range(num_to_avg):
		temp = mylib.noise(steps, (j/2.), fhigh, flow, Gain, h)
		noise = temp
		noise_sum = list(map(add, noise_sum, noise))
		fft_temp = numpy.fft.fft(temp)
		psd_temp = (abs((fft_temp[0:(steps_2)])/steps)**2.)*(1./(fhigh*steps))*2 #temporary one-sided spectrum
		psd_sum = list(map(add, psd_temp, psd_sum)) #element-wise addition
		sys.stdout.write('%s\r' % (((i+1)/(num_to_avg*1.))*100))
		sys.stdout.flush()
	
	#elifs here to divide gammas in to respective arrays
	
	if j == 0:
		psd_avg0 = [n/((num_to_avg*1.)) for n in psd_sum]
		psd_avg0 = [log(p, 10) for p in psd_avg0]
		noise0 = [d/num_to_avg for d in noise_sum]
		print("Finished gamma=0")
		
	elif j == 1:
		psd_avg05 = [n/((num_to_avg*1.)) for n in psd_sum]
		psd_avg05 = [log(p, 10) for p in psd_avg05]
		noise05 = [d/num_to_avg for d in noise_sum]
		print("Finished gamma=0.5")
		
	elif j == 2:
		psd_avg1 = [n/((num_to_avg*1.)) for n in psd_sum]
		psd_avg1 = [log(p, 10) for p in psd_avg1]
		noise1 = [d/num_to_avg for d in noise_sum]
		print("Finished gamma=1")
		
	elif j == 3:
		psd_avg15 = [n/((num_to_avg*1.)) for n in psd_sum]
		psd_avg15 = [log(p, 10) for p in psd_avg15]
		noise15 = [d/num_to_avg for d in noise_sum]
		print("Finished gamma=1.5")
		
	else:
		psd_avg2 = [n/((num_to_avg*1.)) for n in psd_sum]
		psd_avg2 = [log(p, 10) for p in psd_avg2]
		noise2 = [d/num_to_avg for d in noise_sum]
		print("Finished gamma=2.0")
		

print("--- %s seconds ---" % (time.time() - start_time)) #Will give time to compute

Fs = min(fhigh, 0.5/(1./fhigh))                 # Sampling Frequency
#Fs = fhigh
T = 1./Fs                  						# Sampling Period  
L = steps                  						# Length of Signal
t = numpy.arange(steps)*T  						# Time Vector
f = (numpy.arange(L/2.)/L)*Fs

#Linear regression check...

f[0] = 0.001
model = LinearRegression()
test_x = [log(i, 10) for i in f]
test_x = numpy.array(test_x)
model.fit(test_x.reshape((-1,1)), psd_avg0)
print(model.coef_)
model.fit(test_x.reshape((-1,1)), psd_avg05)
print(model.coef_)
model.fit(test_x.reshape((-1,1)), psd_avg1)
print(model.coef_)
model.fit(test_x.reshape((-1,1)), psd_avg15)
print(model.coef_)
model.fit(test_x.reshape((-1,1)), psd_avg2)
print(model.coef_)

# Set up the axes with gridspec
fig = plt.figure(figsize=(15, 10))
grid = plt.GridSpec(10, 6, hspace=0.5, wspace=0.5)

plt.subplot(grid[0:, 2:4])

plt.plot(f, psd_avg0, label=r"$\alpha = 0$")
plt.plot(f, psd_avg05, label=r"$\alpha = 0.5$")
plt.plot(f, psd_avg1, label=r"$\alpha = 1$")
plt.plot(f, psd_avg15, label=r"$\alpha = 1.5$")
plt.plot(f, psd_avg2, label=r"$\alpha = 2$")
plt.title(r"$\dfrac{1}{f^{\alpha}}$ Power Spectral Density", fontsize=16, y=1.03)
plt.legend(loc="upper right", fontsize=14)
pyplot.xscale('log')
plt.xlim(0.05, 500)
#plt.ylim(-140, -40)
plt.xlabel('Hz (log scale)', fontsize=14)
plt.ylabel('Power/Frequency (dB/Hz)', fontsize=14)


plt.subplot(grid[0:2, :2])
plt.title("Noise", fontsize=16, y=1.11)
plt.plot(t, noise0, label=r"$\alpha = 0$", color='#1f77b4')
plt.legend(loc="upper right", fontsize=12)
plt.xticks([])
#plt.yticks([])

plt.subplot(grid[2:4, :2])
plt.plot(t, noise05, label=r"$\alpha = 0.5$", color='#ff7f0e')
plt.legend(loc="upper right", fontsize=12)
plt.xticks([])
#plt.yticks([])

plt.subplot(grid[4:6, :2])
plt.plot(t, noise1, label=r"$\alpha = 1$", color='#2ca02c')
plt.legend(loc="upper right", fontsize=12)
plt.ylabel('Amplitude', fontsize=14)
plt.xticks([])
#plt.yticks([])

plt.subplot(grid[6:8, :2])
plt.plot(t, noise15, label=r"$\alpha = 1.5$", color='#d62728')
plt.legend(loc="upper right", fontsize=12)
plt.xticks([])
#plt.yticks([])

plt.subplot(grid[8:10, :2])
plt.plot(t, noise2, label=r"$\alpha = 2$", color='#9467bd')
plt.legend(loc="upper right", fontsize=12)
plt.xlabel("Time", fontsize=14, labelpad=22)
plt.xticks([])
#plt.yticks([])

plt.show()	
