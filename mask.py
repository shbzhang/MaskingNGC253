import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy import units as u
from scipy.ndimage import convolve
import radio_beam
from spectral_cube import SpectralCube
import cc3d
import tqdm
import matplotlib.pyplot as plt
from reproject import reproject_interp
from astropy.wcs import WCS
from mwispy import DataCube as dc


#c = fits.open('ngc253-h13cn10-cut.fits')[0]
#c = fits.open('ngc253-hcnhp21-cut.fits')[0]
#c = fits.open('ngc253-hcnhp32-cut.fits')[0]
c = fits.open('ngc253-hcnhp43-cut.fits')[0]
#c = fits.open('ngc253-hcnhp54-cut.fits')[0]
#h = c.header
#d = c.data


def getNegRms(data):
	### estimate noise from negative values
	return np.sqrt(np.mean(data[data<0]**2))


def labelCube(data, low, peak):
	### connected components with cc3d
	print('label cc')
	labels = cc3d.connected_components(data>low, connectivity=26)

	print('max=',labels.max())
	mask = np.zeros(labels.shape, dtype=bool)
	for i in tqdm.tqdm(range(1,labels.max()+1)):
		cc = labels==i
		if cc.sum()<19: continue
		if cc.any(axis=(1,2)).sum()<=2: continue
		if data[cc].max()<peak: continue
		mask[cc] = 1
	return mask


def growMask(mask, sigma):
	### grow mask with an ellipsoid
	print('grow mask')
	def TriaxialKernel(ksize, sigma):
		# Create a range from -kernel_size/2 to kernel_size/2
		x, y, z = np.mgrid[-(ksize[0]//2):ksize[0]//2+1,
						   -(ksize[1]//2):ksize[1]//2+1,
						   -(ksize[2]//2):ksize[2]//2+1]
		#kernel = np.exp(-(x**2/sigma[0]**2 + y**2/sigma[1]**2 + z**2/sigma[2]**2) / 2)	    
		#kernel /= kernel.sum()
		kernel = (x**2/sigma[0]**2 + y**2/sigma[1]**2 + z**2/sigma[2]**2) < 1
		return kernel

	ksize = tuple(int(s)*2+3 for s in sigma)
	k = TriaxialKernel(ksize, sigma)
	convMask = convolve(mask, k)
	return (convMask>0).astype(int)


def reprojectMask(mask, target):
	### reproject spatially and spectrally
	print('reproject')
	dcMask = dc.read(mask)
	resMask = dcMask.resample(target.header, largecube=True)

	newdata = np.ndarray(target.data.shape, dtype = float)
	oldwcs = WCS(resMask.header, naxis=2)
	newwcs = WCS(target.header, naxis=2)
	shape_out = target.data.shape[1:]
	for i in range(resMask.shape[0]):
		newdata[i], coverage = reproject_interp((resMask._data[i], oldwcs), newwcs, shape_out = shape_out, order='nearest-neighbor')
	newdata = (newdata>0).astype(int)

	return fits.PrimaryHDU(data=newdata, header=target.header)




### label cube, get mask
if 0:
	c = fits.open('ngc253-hcnhp21-cut.fits')[0]	# use 2-1 as reference as it is not contaminated
	rms = getNegRms(c.data)
	mask = labelCube(c.data, rms*3, rms*5)
	fits.PrimaryHDU(data=mask.astype(int), header=c.header).writeto('mask21.fits', overwrite=True)
else:
	mask = fits.open('mask21.fits')[0].data


### grow mask with kernel
if 0:
	c = fits.open('ngc253-hcnhp21-cut.fits')[0]
	if 1:
		# ordinary mask
		mask = growMask(mask, (3, 17, 17))
		maskHDU = fits.PrimaryHDU(data=mask, header=c.header)
		maskHDU.writeto('mask21_exp925.fits', overwrite=True)
	else:
		# grow a larger mask
		mask = growMask(mask, (9, 25, 25))
		maskHDU = fits.PrimaryHDU(data=mask, header=c.header)
		maskHDU.writeto('mask21_exp925.fits', overwrite=True)


### reproject to align
if 0:
	#use larger mask to include extended emission
	maskHDU = fits.open('mask21_exp925.fits')[0]
	c = fits.open('ngc253-h13cn10-cut.fits')[0]
	rpjHDU = reprojectMask(maskHDU, c)
	rpjHDU.writeto('mask21_exp925_rpj10.fits', overwrite=True)

	#use ordinary mask for more compact lines
	maskHDU = fits.open('mask21_exp317.fits')[0]
	c = fits.open('ngc253-hcnhp32-cut.fits')[0]
	rpjHDU = reprojectMask(maskHDU, c)
	rpjHDU.writeto('mask21_exp317_rpj32.fits', overwrite=True)

	c = fits.open('ngc253-hcnhp43-cut.fits')[0]
	rpjHDU = reprojectMask(maskHDU, c)
	rpjHDU.writeto('mask21_exp317_rpj43.fits', overwrite=True)

	c = fits.open('ngc253-hcnhp54-cut.fits')[0]
	rpjHDU = reprojectMask(maskHDU, c)
	rpjHDU.writeto('mask21_exp317_rpj54.fits', overwrite=True)



if 0:
	# remove wrong GMC 3,5 in 4-3, it is from other line
	c = fits.open('ngc253-hcnhp43-cut.fits')[0]
	rms = getNegRms(c.data)
	labels = cc3d.connected_components(c.data>rms*14, connectivity=26)
	# find the label of wrong GMC
	print('Wrong label: ', labels[83, 235, 362])
	empty = ((labels==4) | (labels==3)).astype(int)

	def Circle3DKernel(ksize, sigma):
		# Create a range from -kernel_size/2 to kernel_size/2
		x, y, z = np.mgrid[-(ksize[0]//2):ksize[0]//2+1,
						   -(ksize[1]//2):ksize[1]//2+1,
						   -(ksize[2]//2):ksize[2]//2+1]
		kernel = (x**2/sigma[0]**2 + y**2/sigma[1]**2 + z**2/sigma[2]**2) < 1
		return kernel

	sigma = [1.1, 7, 7]
	ksize = tuple(int(s)*2+3 for s in sigma)
	k = Circle3DKernel(ksize, sigma)
	emptyMask = convolve(empty, k)

	hdu = fits.open('mask21_exp317_rpj43.fits')[0]
	print(hdu.data.shape, emptyMask.shape)
	hdu.data[emptyMask.astype(bool)] = 0
	hdu.writeto('mask21_exp317_rpj43_remove.fits', overwrite=True)
	#fits.PrimaryHDU(data=mask.astype(int), header=c.header).writeto('mask43.fits', overwrite=True)

	c = fits.open('ngc253-hcnhp43-cut.fits')[0]
	hdu = fits.open('mask21_exp317_rpj43_remove.fits')[0]
	plt.imshow(-np.sum(hdu.data * c.data, axis=0)*hdu.header['CDELT3'], origin='lower', cmap='rainbow', vmin=-0.7, vmax=2)
	#plt.plot([336, 320], [257,270], 'r+')
	#plt.show()


### quicklook results
if 1:
	fig, ax = plt.subplots(nrows=2, ncols=3)#, sharex=True, sharey=True)
	ax = ax.flatten()

	c = fits.open('ngc253-h13cn10-cut.fits')[0].data
	m = fits.open('mask10.fits')[0].data
	ma = fits.open('mask21_exp925_rpj10.fits')[0].data
	np.save('mask10.npy', ma.astype(bool))
	ax[0].imshow(np.nansum(c*ma, axis=0), origin='lower', aspect='equal', cmap='rainbow')
	ax[0].contour(np.any(ma, axis=0), levels=[0.5], colors='r')

	c = fits.open('ngc253-hcnhp21-cut.fits')[0].data
	m = fits.open('mask21.fits')[0].data
	ma = fits.open('mask21_exp317.fits')[0].data
	np.save('mask21.npy', ma.astype(bool))
	ax[1].imshow(np.nansum(c*ma, axis=0), origin='lower', aspect='equal', cmap='rainbow')
	ax[1].contour(np.any(ma, axis=0), levels=[0.5], colors='r')

	c = fits.open('ngc253-hcnhp32-cut.fits')[0].data
	m = fits.open('mask32.fits')[0].data
	ma = fits.open('mask21_exp317_rpj32.fits')[0].data
	np.save('mask32.npy', ma.astype(bool))
	ax[2].imshow(np.nansum(c*ma, axis=0), origin='lower', aspect='equal', cmap='rainbow')
	ax[2].contour(np.any(ma, axis=0), levels=[0.5], colors='r')

	c = fits.open('ngc253-hcnhp43-cut.fits')[0].data
	m = fits.open('mask43.fits')[0].data
	ma = fits.open('mask21_exp317_rpj43_remove.fits')[0].data
	np.save('mask43.npy', ma.astype(bool))
	ax[3].imshow(np.nansum(c*ma, axis=0), origin='lower', aspect='equal', cmap='rainbow')
	ax[3].contour(np.any(ma, axis=0), levels=[0.5], colors='r')

	c = fits.open('ngc253-hcnhp54-cut.fits')[0].data
	m = fits.open('mask54.fits')[0].data
	ma = fits.open('mask21_exp317_rpj54.fits')[0].data
	np.save('mask54.npy', ma.astype(bool))
	ax[4].imshow(np.nansum(c*ma, axis=0), origin='lower', aspect='equal', cmap='rainbow')
	ax[4].contour(np.any(ma, axis=0), levels=[0.5], colors='r')

plt.show()