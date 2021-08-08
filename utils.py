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

	# Create copy of input angle for output name
	phiIn = phi0

	# Offset in angle to fix the phi such that white brushes line up with angle
	phi0 = (phi0 + 113)%180

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
		whole = str(phiIn).split('.')[0].zfill(3)

		# Decimal part of angle
		decimal = str(phiIn).split('.')[-1]

		# Create name for output image
		maskPath = 'masks/r_' + str(ra) + '_a_' + whole + '_' + decimal + '.png'

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
def find_angle(defectZoom, ra):

	# Where masks are stored, check for those of the correct radius
	maskPaths = glob.glob('masks/r_' + str(ra) + '_a_*.png')

	# Storing totals
	tots = []

	# Storing angles
	angles = []

	# Go through every mask
	for i, maskPath in enumerate(maskPaths):

		# String describing angle
		anStr = maskPath.split('\\')[-1].split('.')[0][-5:]

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

	if True:

		# Showing defect and selected template
		plt.imshow(defectZoom, cmap='gray')
		plt.show()
		plt.imshow(Image.open(maskPaths[index]), cmap='gray')
		plt.show()

		# Showing plot of totals vs angles
		plt.plot(angles, tots)
		plt.scatter(bestAngle, np.amax(tots), color='r')
		plt.show()

	# Return result
	return bestAngle



# Add a border to a numpy array
def add_border(array, width, value=0):

	# Matrix to store output
	mat = []

	# Create every row
	for i in range(len(array) + 2*width):

		# List to store values in a row
		row = []

		# Go through every column in row
		for j in range(len(array[0]) + 2*width):

			# If it exceeds the dimensions of the original
			if i - width < 0 or j - width < 0 or i - width >= len(array) or j - width >= len(array[0]):

				# Add the value specified
				row.append(value)

			# Otherwise copy over the image
			else:

				# Image pixel value
				row.append(array[i-width][j-width])

		# Append the row as numpy array to final matrix
		mat.append(np.asarray(row))

	# Return the matrix as numpy array
	return np.asarray(mat)



## Find defect angles in an image
def image_angles(imPath):

	print('\nDetecting defect angles in ' + imPath)

	# Get path to label file corresponding to image
	labelPath = imPath.replace('images', 'labels').replace('tif', 'txt')

	# If the corresponding label file exists,
	if len(glob.glob(labelPath)) == 1:

		# Create path to output label file
		outPath = labelPath.replace('input', 'output')

		# Read the image as numpy array
		image = np.asarray(Image.open(imPath))

		# Read labels as numpy array
		labelMat = np.loadtxt(labelPath)

		# Output matrix
		outMat = []

		# Pad numpy image with extra 0 values
		imageFat = add_border(image, 100, 0)

		# Store whether it is 0 or something else
		zeros = 0
		nonzeros = 0

		# For every angle in the defect
		for i, defect in enumerate(labelMat):

			# If the length of the 'defect' label is nonzero
			if not np.isscalar(defect) and len(defect) > 0:

				# Find center of detected defect
				xc = int(defect[1]*image.shape[1])
				yc = int(defect[2]*image.shape[0])

				## Now find the nearest neighbour to this defect
				# If at least two defects exist
				if len(labelMat) >= 2:

					# Create a temporary minDist2
					minDist2 = 1e20

					# Go through every defect
					for j, defectPair in enumerate(labelMat):

						# If the defect is not the same as the pair
						if i != j:

							# Also, if the length of the acquired defect is correct
							if not np.isscalar(defectPair) and len(defectPair) > 0:

								# Find center of detected defect pair
								xcp = int(defectPair[1]*image.shape[1])
								ycp = int(defectPair[2]*image.shape[0])

								# Find the distance squared
								d2 = (xcp - xc)**2 + (ycp - yc)**2

								# If this is lesser than the last minimum
								if d2 < minDist2:

									# Set this as the new minimum distance squared
									minDist2 = d2

								# Clear variables
								xcp, ycp, d2 = None, None, None

				# Acquired radius of mask template
				ra = int(np.sqrt(minDist2)/4) + 2

				# If the radius is too big, bring it down to 100
				if ra > 100:
					ra = 100

				# Slice out the defect (including pad)
				defArray = imageFat[yc - ra + 100:yc + ra + 101, xc - ra + 100:xc + ra + 101]

				# Get mean value
				mean = np.mean(defArray)

				# Increase contrast
				defArray = (defArray - mean)*100 + mean

				# Get the best defect angle
				bestAngle = find_angle(defArray, ra)

				# Store as zero if zero, otherwise if not
				if int(bestAngle) <= 90:
					zeros = zeros + 1
				else:
					nonzeros = nonzeros + 1

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

		# Display how many were zeros
		print('Under 90: ' + str(zeros) + ', Over 90: ' + str(nonzeros))
		if (zeros+nonzeros) > 0:
			print('Ratio (90/all): ' + str(zeros/(zeros+nonzeros)))

		# Pause if it is too big
		if (zeros+nonzeros) > 0 and zeros/(zeros+nonzeros) > 0.5:
			input('Very big')



if __name__ == '__main__':
	
	maskTemplate(0, ra=10, display=True)