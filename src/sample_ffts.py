#The purpose of this program is to quickly evaluate and test different DFT routines in Python for comparison to the slightly more *ahem* arcane implementation of FFTW

import scipy
import math
import numpy as np
import matplotlib.pyplot as pyplot

size = 20  #Array size for functions


def main():
    f = []
    print np.arange(size)
    for i in np.arange(size):
        #f.append(5)                             #Contant function
        #f.append(math.sin(2*math.pi*i/size))    #Single-frequency sine wave
        f.append(math.sin(2*math.pi*i/size) + math.sin(10*math.pi*i/size))   #Multiple sine waves
    #pyplot.plot(2*math.pi*np.arange(size)/size, f)
    pyplot.plot(np.arange(size), f)
    pyplot.show()
    
    npf = np.array(f)
    print npf
    npf_fft = np.fft.fft(npf)
    print npf_fft
    #pyplot.plot(2*math.pi*np.arange(size)/size, np.imag(npf_fft), 'b')
    #pyplot.plot(2*math.pi*np.arange(size)/size, np.real(npf_fft), 'r')
    #pyplot.plot(2*math.pi*np.arange(size)/size, np.abs(npf_fft), 'k')
    pyplot.plot(np.arange(size), np.imag(npf_fft), 'b')
    pyplot.plot(np.arange(size), np.real(npf_fft), 'r')
    pyplot.plot(np.arange(size), np.abs(npf_fft), 'k')
    pyplot.show()
    
    npf_fft_ifft = np.fft.ifft(npf_fft)
    print npf_fft_ifft
    #pyplot.plot(2*math.pi*np.arange(size)/size, np.real(npf), 'b')
    #pyplot.plot(2*math.pi*np.arange(size)/size, np.real(npf_fft_ifft), 'r')
    pyplot.plot(np.arange(size), np.real(npf), 'b')
    pyplot.plot(np.arange(size), np.real(npf_fft_ifft), 'r')
    pyplot.show()
    
if __name__ == '__main__':
    main()