### Importing packages
import cv2
import numpy as np
import matplotlib.pyplot as plt
import trackpy as tp
import glob
import shutil
import os

from scipy.signal import fftconvolve

### Importing useful scripts
from utils import *

print('Imported packages.')

def normxcorr2(template, image, mode="full"):
	"""
	Input arrays should be floating point numbers.
	:param template: N-D array, of template or filter you are using for cross-correlation.
	Must be less or equal dimensions to image.
	Length of each dimension must be less than length of image.
	:param image: N-D array
	:param mode: Options, "full", "valid", "same"
	full (Default): The output of fftconvolve is the full discrete linear convolution of the inputs. 
	Output size will be image size + 1/2 template size in each dimension.
	valid: The output consists only of those elements that do not rely on the zero-padding.
	same: The output is the same size as image, centered with respect to the ‘full’ output.
	:return: N-D array of same dimensions as image. Size depends on mode parameter.
	"""

	# If this happens, it is probably a mistake
	if np.ndim(template) > np.ndim(image) or \
			len([i for i in range(np.ndim(template)) if template.shape[i] > image.shape[i]]) > 0:
		print("normxcorr2: TEMPLATE larger than IMG. Arguments may be swapped.")

	template = template - np.mean(template)
	image = image - np.mean(image)

	a1 = np.ones(template.shape)
	# Faster to flip up down and left right then use fftconvolve instead of scipy's correlate
	ar = np.flipud(np.fliplr(template))
	out = fftconvolve(image, ar.conj(), mode=mode)
	
	image = fftconvolve(np.square(image), a1, mode=mode) - \
			np.square(fftconvolve(image, a1, mode=mode)) / (np.prod(template.shape))

	# Remove small machine precision errors after subtraction
	image[np.where(image < 0)] = 0

	template = np.sum(np.square(template))
	out = out / np.sqrt(image * template)

	# Remove any divisions by 0 or very close to 0
	out[np.where(np.logical_not(np.isfinite(out)))] = 0
	
	return out



def create_cc_masks():

	# 'Radius' of the defect ((mask dimension - 1)/2)
	ra = 20

	# Number of brushes in defect template
	# (4 for crossed polarizers, 2 for decrossed)
	nBrushes = 2

	# Minimum template orientation angle (degrees)
	phiMin = 0

	# Maximum template orientation angle (degrees)
	# (90 for crossed polarizers, 180 for decrossed)
	phiMax = 135

	# Number of masks to generate
	nMasks = 16

	# Display the final result?
	dispResults = False

	# Display the individual components?
	dispComp = False

	# Create mask output directory
	if len(glob.glob('masks/')) != 0:
		shutil.rmtree('masks/')
	os.mkdir('masks/')

	# Store various angles in the range
	tAngle = np.linspace(phiMin, phiMax, num=nMasks)

	# Create a mask for each 
	for i, phi0 in enumerate(tAngle):

		# Generate a mask
		maskTemplate(phi0, ra=ra, nBrushes=nBrushes, dispComp=dispComp, dispResults=dispResults)

def main():

	print('\nRunning main.')

	# imPath = 'input/' + '2021-06-20_17-52-22-t_17.tiff'
	imPath = 'input/' + 'r2_0110.tif'

	image = cv2.imread(imPath, cv2.IMREAD_GRAYSCALE)/255
	template = cv2.imread('kernel.png', cv2.IMREAD_GRAYSCALE)/255

	print('\nPerforming correlate.')
	corr = normxcorr2(template, image)

	print('\nThresholding result.')
	thresh = (corr > 0.4) * corr

	print('\nLocating bright features')
	f = tp.locate(thresh, 5, invert=False)

	x = f['x'] - template.shape[0]//2
	y = f['y'] - template.shape[1]//2

	fig = plt.figure(figsize=(8, 8))

	fig.add_subplot(1, 3, 1)
	plt.scatter(x, y, color='r', s=2)
	plt.imshow(image, cmap='gray')

	fig.add_subplot(1, 3, 2)
	plt.imshow(corr, cmap='gray')

	# Thresholded corr
	fig.add_subplot(1, 3, 3)
	plt.imshow(thresh, cmap='gray')

	plt.show()

if __name__ == '__main__':

	# main()
	create_cc_masks()