### Importing packages
import cv2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as ps
import trackpy as tp
import glob
import shutil
import os
from PIL import Image
import shutil

from scipy.signal import fftconvolve

### Importing useful scripts
from utils import *

print('Imported packages.')

def create_cc_masks():

	# Number of brushes in defect template
	# (4 for crossed polarizers, 2 for decrossed)
	nBrushes = 2

	# Minimum template orientation angle (degrees)
	phiMin = 0

	# Maximum template orientation angle (degrees)
	# (90 for crossed polarizers, 180 for decrossed)
	phiMax = 180

	# Number of masks to generate
	nMasks = 181

	# Display the final result?
	dispResults = False

	# Display the individual components?
	dispComp = False

	# Save the results?
	saveMask = True

	# Recreate mask output only if results are to be saved
	if saveMask:

		# Create mask output directory
		if len(glob.glob('masks/')) != 0:
			shutil.rmtree('masks/')
		os.mkdir('masks/')

	# Store various angles in the range
	tAngle = np.linspace(phiMin, phiMax, num=nMasks)

	# Create a mask of every radius between 1 and 100
	for ra in range(1, 101):

		print(str(ra)+'%')

		# Create a mask for each angle
		for i, phi0 in enumerate(tAngle):

			# Generate a mask
			maskTemplate(phi0, ra=ra, nBrushes=nBrushes, dispComp=dispComp, dispResults=dispResults, saveMask=saveMask)



## Main functioning of script
def main():

	# If output folder exists
	if len(glob.glob('output/')) > 0:
		shutil.rmtree('output/')

	# Create output folder
	os.mkdir('output/')

	# Create output labels
	os.mkdir('output/labels/')

	# Read input images
	imagePaths = glob.glob('input/images/*.tif')

	# Perform angle extraction on each image
	for i, imagePath in enumerate(imagePaths):

		imagePath = 'input/images/r2_4000.tif'

		# Extract angles
		image_angles(imagePath)



if __name__ == '__main__':

	# create_cc_masks()

	main()