## Importing useful packages
import numpy as np
import glob
from PIL import Image
import matplotlib.pyplot as plt


## Cross corellate image with a template
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
	
	# Return cross-correlation result
	return out



## Take in the rotation angle of a defect and create defect mask array
def maskTemplate(phi0, ra=5, nBrushes=2, dispComp=False, dispResults=False, saveMask=True):

	# Conversion factor degrees to radians
	dtoR = np.pi/180

	# Create template x+iy and modify it to  ((x+iy)/(SQRT(x^2 + y^2)))^(-4) for 4-brush pattern
	txPixels = 2*ra + 1
	tyPixels = 2*ra + 1

	# Complex number arrays
	real_z = np.zeros((txPixels, tyPixels))
	imag_z = np.zeros((txPixels, tyPixels))

	# Indexed vector [-ra, ... -1, 0, 1, ... ra]
	vRange = np.linspace(-ra, ra, num=2*ra+1)

	# Replace row/column in complex number arrays
	for i in range(2*ra + 1):

		# Fill rows with this vector ("x")
		real_z[:, i] = vRange

		# Fill columns with this vector ("y")
		imag_z[i, :] = vRange

	# Force center values to unity
	real_z[ra, ra] = 1
	imag_z[ra, ra] = 1

	# Gives x + iy
	z0  = real_z + imag_z*1j

	# x + iy rotated  by -90 degrees (y - ix)   (was +Re, -Im but this rotates the patern clockwise ...)
	z90 = imag_z - real_z*1j

	# Combine the orthogonal complex fields at angle phi0
	z = z0*np.cos(phi0*dtoR) + z90*np.sin(phi0*dtoR)

	# Normalize  to unit vectors so have exp[i(phi - phi0)]
	template = z/np.sqrt(real_z**2 + imag_z**2)

	# Invert to create a sharp 4-fold cross
	# template = 1./template^4

	# Invert to create a sharp cross with nbrushes arms
	# jem april 2009
	template = 1/template**nBrushes

	# Zero the elements outside radius ra
	template = (np.sqrt(real_z**2 + imag_z**2) < ra)*template

	# Compute real and imaginary parts of the template
	# (although seems none of these was passed back to the main program in the original code - jem)
	realT = (np.real(template) + 1)
	imaginaryT = (np.imag(template) + 1)

	# Direct sum of the real/imaginary parts (mean)
	bothT = (realT + imaginaryT)

	# Normalize bothT
	bothT = ((bothT/np.amax(bothT))*255).astype(int)

	# If the user wanted to display each parts
	if dispComp:

		# Displaying the real part

		# Clear plot
		plt.clf()

		# Plot real part as image
		plt.imshow(realT, cmap='gray')

		# Add colorbar
		plt.colorbar()

		# Show plot
		plt.show()

		# Displaying the imaginary part

		# Clear plot
		plt.clf()

		# Plot real part as image
		plt.imshow(imaginaryT, cmap='gray')

		# Add colorbar
		plt.colorbar()

		# Show plot
		plt.show()

	# If the user wanted to display final results
	if dispResults:

		# Displaying resultant plot

		# Clear plot
		plt.clf()

		# Plot real part as image
		plt.imshow(bothT, cmap='gray')

		# Add colorbar
		plt.colorbar()

		# Show plot
		plt.show()

	# If user wants to save the masks
	if saveMask:

		# Whole part of angle
		whole = str(phi0).split('.')[0].zfill(3)

		# Decimal part of angle
		decimal = str(phi0).split('.')[-1]

		# Create name for output image
		maskPath = 'masks/' + whole + '_' + decimal + '.png'

		# Numpy array to image
		maskIm = Image.fromarray(np.uint8(bothT))

		# Save mask
		maskIm.save(maskPath)



## Elementwise sum for two arrays
def dot_prod(array1, array2):

	# Multiply each element
	prod = np.multiply(array1, array2)

	# Get total
	total = np.sum(prod)

	# Return result
	return total



## Find the best angle the defect correlates to
def find_angle(defectZoom):

	# Where masks are stored
	maskPaths = glob.glob('masks/*.png')

	# Storing totals
	tots = []

	# Storing angles
	angles = []

	# Go through every mask
	for i, maskPath in enumerate(maskPaths):

		# String describing angle
		anStr = maskPath.split('\\')[-1].split('.')[0]

		# Find angle and add to array
		angles.append(float(anStr.replace('_', '.')))

		# Read mask template
		template = np.asarray(Image.open(maskPath))

		# Find elementwise sum
		tots.append(dot_prod(defectZoom, template))

	# Convert totals list to numpy array
	tots = np.asarray(tots)

	# Index of largest sum
	index = np.where(tots == np.amax(tots))[0][0]

	# Get best angle from list
	bestAngle = angles[index]

	# # Showing defect and selected template
	# plt.imshow(defectZoom, cmap='gray')
	# plt.show()
	# plt.imshow(Image.open(maskPaths[index]), cmap='gray')
	# plt.show()

	# # Showing plot of totals vs angles
	# plt.plot(angles, tots)
	# plt.scatter(bestAngle, np.amax(tots), color='r')
	# plt.show()

	# Return result
	return bestAngle



## Find defect angles in an image
def image_angles(imPath):

	print('\nDetecting defect angles in ' + imPath)

	# Get path to label file corresponding to image
	labelPath = imPath.replace('images', 'labels').replace('tif', 'txt')

	# Create path to output label file
	outPath = labelPath.replace('input', 'output')

	# Read the image as numpy array
	image = np.asarray(Image.open(imPath))

	# Read labels as numpy array
	labelMat = np.loadtxt(labelPath)

	# Output matrix
	outMat = []

	# Get semi width and length of mask
	testMask = np.asarray(Image.open(glob.glob('masks/*.png')[0]))
	sw = (len(testMask) - 1)//2

	# For every angle in the defect
	for i, defect in enumerate(labelMat):

		# Find center of detected defect
		xc = int(defect[1]*len(image[0]))
		yc = int(defect[2]*len(image))

		# Slice out the defect
		defArray = image[yc - sw:yc + sw+1, xc - sw:xc + sw+1]

		# Get mean value
		mean = np.mean(defArray)

		# Increase contrast
		defArray = (defArray - mean)*100 + mean

		# Get the best defect angle
		bestAngle = find_angle(defArray)

		# Output matrix line
		line = np.copy(defect)

		# Now add the angle
		line = np.append(line, bestAngle)

		# Add this to the output array
		outMat.append(line)

	# Convert outarray to numpy array
	outMat = np.asarray(outMat)

	# Save output matrix as txt
	np.savetxt(outPath, outMat)



if __name__ == '__main__':
	
	maskTemplate(0, ra=10, display=True)