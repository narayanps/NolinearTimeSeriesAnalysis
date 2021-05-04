#AUTHOR : NARAYAN P SUBRAMANIYAM
import numpy


def UnivariateSurrogates(data_f,MaxIter):
    
    xs=data_f.copy()
    xs.sort() #sorted amplitude stored
    pwx=numpy.abs(numpy.fft.fft(data_f)) # amplitude of fourier transform of orig
    
    data_f.shape = (-1,1)
    xsur = numpy.random.permutation(data_f) #random permutation as starting point
    xsur.shape = (1,-1)
    
    for i in range(MaxIter):
        fftsurx = pwx*numpy.exp(1j*numpy.angle(numpy.fft.fft(xsur)))
        xoutb = numpy.real(numpy.fft.ifft(fftsurx))
        ranks = xoutb.argsort(axis=1)
        xsur[:,ranks] = xs
    return(xsur)


def BivariateSurrogates(data_f,MaxIter):
    
     M = data_f.shape[1]

     sorted_original = data_f.copy()
     sorted_original.sort(axis=1)

     fourier_transform = numpy.fft.fft(data_f)
     original_fourier_amps = numpy.abs(fourier_transform)
     original_fourier_phase = numpy.angle(fourier_transform)
     phix = original_fourier_phase
     niterations = 120
     xsur = numpy.random.permutation(data_f).transpose()
     
     for i in xrange(niterations):
          phisx = numpy.angle(numpy.fft.fft(xsur))
          alpha = numpy.zeros((M,1))
          for cc in range(1,int(numpy.ceil(M/2)),1):
               alpha[cc] = numpy.arctan(numpy.sum(numpy.sin(phisx[:,cc]-phix[:,cc]))/numpy.sum(numpy.cos(phisx[:,cc]-phix[:,cc])))
               if numpy.sum(numpy.cos(alpha[cc]+phix[:,cc]-phisx[:,cc])) < 0:
                   alpha[cc] = alpha[cc]-numpy.pi
               alpha[M-1-cc+1] = -alpha[cc]
          alpha.shape = (1,-1)
          tmpx = numpy.tile(alpha,(2,1))
          phases = phix + tmpx
          fftsurx = original_fourier_amps*numpy.exp(1j * phases)
          xoutb = numpy.real(numpy.fft.ifft(fftsurx,axis=1))
          ranks = xoutb.argsort(axis=1)
          for k in xrange(sorted_original.shape[0]):
              xsur[k,ranks[k,:]] = sorted_original[k,:]         
     return(xsur)
     
     
def UnivariateSurrogatesTFT(data_f,MaxIter,fc):
    xs=data_f.copy()
    xs.sort() #sorted amplitude stored
    pwx=numpy.abs(numpy.fft.fft(data_f)) # amplitude of fourier transform of orig
    phi = numpy.angle(numpy.fft.fft(data_f))
    Len=phi.shape[1]
    data_f.shape=(-1,1)
    xsur = numpy.random.permutation(data_f) #random permutation as starting point
    xsur.shape = (1,-1)
    
    Fc =  numpy.round(fc*data_f.shape[0])
    for i in range(MaxIter):
        phi_surr = numpy.angle(numpy.fft.fft(xsur))
        phi_surr[0,1:Fc] = phi[0,1:Fc]
        phi_surr[0,Len-Fc+1:Len] = phi[0,Len-Fc+1:Len]
        phi_surr[0,0] = 0.0
        phi_surr[0,Len/2] = 0.0 
              
        fftsurx = pwx*numpy.exp(1j*phi_surr)
        xoutb = numpy.real(numpy.fft.ifft(fftsurx))
        ranks = xoutb.argsort(axis=1)
        xsur[:,ranks] = xs
    return(xsur)

def BivariateSurrogatesTFT(data_f,MaxIter,fc):
    
     M = data_f.shape[1]
     Fc =  numpy.round(fc*M)
     sorted_original = data_f.copy()
     sorted_original.sort(axis=1)

     fourier_transform = numpy.fft.fft(data_f)
     original_fourier_amps = numpy.abs(fourier_transform)
     original_fourier_phase = numpy.angle(fourier_transform)
     phix = original_fourier_phase
     xsur = numpy.random.permutation(data_f.transpose()).transpose()
     for i in xrange(MaxIter):
          phisx = numpy.angle(numpy.fft.fft(xsur))
          
          phisx[:,1:Fc] = phix[:,1:Fc]
          phisx[:,M-Fc+1:M] = phix[:,M-Fc+1:M]
          phisx[0,0] = 0
          phisx[1,0] = 0
          phisx[:,M/2] = phix[:,M/2]

          
          alpha = numpy.zeros((M,1))
          for cc in range(1,int(numpy.ceil(M/2)),1):
               alpha[cc] = numpy.arctan(numpy.sum(numpy.sin(phisx[:,cc]-phix[:,cc]))/numpy.sum(numpy.cos(phisx[:,cc]-phix[:,cc])))
               if numpy.sum(numpy.cos(alpha[cc]+phix[:,cc]-phisx[:,cc])) < 0:
                   alpha[cc] = alpha[cc]-numpy.pi
               alpha[M-1-cc+1] = -alpha[cc]
          alpha.shape = (1,-1)
          tmpx = numpy.tile(alpha,(2,1))
          phases = phix + tmpx
          fftsurx = original_fourier_amps*numpy.exp(1j * phases)
          xoutb = numpy.real(numpy.fft.ifft(fftsurx,axis=1))
          ranks = xoutb.argsort(axis=1)
          for k in xrange(sorted_original.shape[0]):
              xsur[k,ranks[k,:]] = sorted_original[k,:]         
     return(xsur)
