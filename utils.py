## Importing useful packages
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt



## Take in the rotation angle of a defect and create defect mask array
def maskTemplate(phi0, ra=5, nBrushes=2, dispComp=False, dispResults=False):

	# Conversion factor degrees to radians
	dtoR = np.pi/180

	# create template x+iy and modify it to  ((x+iy)/(SQRT(x^2 + y^2)))^(-4) for 4-brush patter
	txPixels = 2*ra + 1
	tyPixels = 2*ra + 1

	# Complex number arrays
	real_z = np.zeros((txPixels, tyPixels))
	imag_z = np.zeros((txPixels, tyPixels))

	# indexed vector [-ra, ... -1, 0, 1, ... ra]
	vRange = np.linspace(-ra, ra, num=2*ra+1)

	# Replace row/column in complex number arrays
	for i in range(2*ra + 1):

		# fill rows with this vector ("x")
		real_z[:, i] = vRange

		# fill columns with this vector ("y")
		imag_z[i, :] = vRange

	# force center values to unity
	real_z[ra, ra] = 1
	imag_z[ra, ra] = 1

	# gives x + iy
	z0  = real_z + imag_z*1j

	# x+iy rotated  by -90 degrees (~ y -ix)   (was +Re, -Im but this rotates the patern clockwise ...)
	z90 = imag_z - real_z*1j

	# combine the orthogonal complex fields at angle phi0
	z = z0*np.cos(phi0*dtoR) + z90*np.sin(phi0*dtoR)

	# normalize  to unit vectors so have exp[i(phi - phi0)]
	template = z/np.sqrt(real_z**2 + imag_z**2)

	# invert to create a sharp 4-fold cross
	# template = 1./template^4

	# invert to create a sharp cross with nbrushes arms
	# jem april 2009
	template = 1/template**nBrushes

	# zero the elements outside radius ra
	template = (np.sqrt(real_z**2 + imag_z**2) < ra)*template

	# compute real and imaginary parts of the template
	# (although seems none of these was passed back to the main program in the original code - jem)
	realT = ((np.real(template) + 1)/2)*255
	imaginaryT = ((np.imag(template) + 1)/2)*255

	# Direct sum of the real/imaginary parts (mean)
	bothT = (realT + imaginaryT)/2

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

	# Create name for output image
	maskPath = 'masks/' + str(phi0).replace('.', '_') + '.png'

	# Numpy array to image
	maskIm = Image.fromarray(np.uint8(bothT))

	# Save mask
	maskIm.save(maskPath)


if __name__ == '__main__':
	
	maskTemplate(0, ra=10, display=True)