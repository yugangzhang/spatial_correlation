from scipy.ndimage.filters import gaussian_filter
import numpy as np
# not great to do but it makes debugging life so much easier
from matplotlib.pyplot import * 
ion()

def circavg(img,x0=None,y0=None,mask=None,SIMG=None):
    ''' Compute a circular average of the data. 
        If SIMG is not null, put the data into this image.
    '''
    if(mask is None):
        mask = np.ones(img.shape)
    dimy,dimx = img.shape
    if(x0 is None):
        x0 = dimx/2
    if(y0 is None):
        y0 = dimy/2
    img = img.ravel()

    pixellist = np.where(mask.ravel()==1)

    x = np.arange(dimx) - x0
    y = np.arange(dimy) - y0
    X,Y = np.meshgrid(x,y)
    R = np.sqrt(X**2 + Y**2).ravel()
    Rd = (R+.5).astype(int).ravel()
    
    noperR = np.bincount(Rd.ravel()[pixellist]).astype(float)
    w = np.where(noperR != 0)

    Rvals = np.bincount(Rd.ravel()[pixellist],weights=R.ravel()[pixellist])[w]/noperR[w]
    Ivals = np.bincount(Rd.ravel()[pixellist],weights=img.ravel()[pixellist])[w]/noperR[w]
    if(SIMG is not None):
        np.copyto(SIMG.ravel(),np.interp(R,Rvals,Ivals))
    return Rvals, Ivals

def mkbox(pxlst,dims):
    ''' Make the image box with mask from the pixels selected and the known imgwidth'''

    imgwidth = dims[0]
    #x and y are flipped???
    #matrix notation!!!
    pixely = pxlst%imgwidth
    pixelx = pxlst//imgwidth

    minpixelx = np.min(pixelx)
    minpixely = np.min(pixely)
    maxpixelx = np.max(pixelx)
    maxpixely = np.max(pixely)

    widthx = maxpixelx - minpixelx + 1
    widthy = maxpixely - minpixely + 1

    oldimg = np.zeros(dims)
    oldimg.ravel()[pxlst] = 1
    mask = np.copy(oldimg[minpixelx:maxpixelx+1,minpixely:maxpixely+1])
    subpxlst = np.where(mask.ravel() == 1)
    print("min/maxpixelx: {},{};min/maxpixely: {},{}".format(minpixelx,maxpixelx,minpixely,maxpixely))

    return subpxlst,mask
    
def autocorr(img1,img2,mask=None):
    '''Compute the autocorrelation of two images.
        Right now does not take mask into account.
        todo: take mask into account (requires extra calculations)
    '''
    if(mask is not None):
        img1 *= mask
        img2 *= mask
    imgc = fftshift(ifft2(fft2(img1)*np.conj(fft2(img2))).real)
    #imgc /= (img1.shape[0]*img1.shape[1])**2
    if(mask is not None):
        maskc = autocorr(mask,mask)
        w = np.where(maskc > 0)[0]
        if(len(w) > 0):
            imgc[w] /= maskc[w] 
        w = np.where(maskc == 0)
        if(len(w)):
            imgc[w] *= 0
    return imgc


# my own very crude simulation of a film with a certain correlation length
# i.e. kapton film (long range variations, leads to Porod's law)
# meant to test acorr since we don't have good data yet

# coordinate system, change to Q in pixels (doesn't matter)
N = 500
x = np.linspace(-N/2,N/2,N)
y = np.linspace(-N/2,N/2,N)

X,Y = np.meshgrid(x,y)

R1 = np.sqrt(X**2 + Y**2)
Q = np.copy(R1)
PHI = np.arctan2(Y,X)

# make a mask that selects a q wedge
mask = np.zeros((N,N))
w = np.where((Q > 90)*(Q < 100)*(PHI > .1)*(PHI < .4))
mask[w] = 1

x0,y0 = 10,20

# the illuminating Gaussian beam
GAUSSBEAM = np.exp(-R1**2/2./100**2)

# random sample
img1 = np.random.random((N,N))
# but it has a correlation length
img1 = gaussian_filter(img1,1)
# illuminate it with a beam
img1 *= GAUSSBEAM
# Fourier transform (won't be SAXS at high q with Ewald sphere)
simg1 = np.abs(fftshift(fft2(img1)))**2

# now azimuthally average to obtain the average sturcture factor
SIMG = np.zeros(img1.shape)
sqx,sqy = circavg(simg1,SIMG=SIMG)

#simg1 /= SIMG

# this piece makes the subimages and masks, if your sub selection is a
# square/rectangle, then you just make a square mask with ones everywhere
# however, from my experience, wedges usually yield better results
# so subpxlst can also be square indexing. don't use ravel and just index:
# subimg1 = simg1[40:100,40:100] for example etc
# just make a mask of 1s so for ex: mask = subimg1*0 + 1

dims = SIMG.shape
pxlst = np.where(mask.ravel())[0]
subpxlst,submask = mkbox(pxlst,dims)
subimg1 = np.zeros(submask.shape)
subimg1.ravel()[subpxlst] = simg1.ravel()[pxlst]

subimg1 /= np.average(subimg1.ravel()[subpxlst])

# the autocorrelation
imgc = autocorr(subimg1,subimg1,mask=submask)


figure(0);cla();
imshow(simg1*mask)
#clim(0,1)

figure(1);cla()
imshow(imgc)

figure(2);cla()
plot(imgc[15])
plot(imgc[:,8])
