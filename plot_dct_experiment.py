import matplotlib.pyplot as plt
import numpy as np

import math
from scipy.fftpack import idct, dct
from mpl_toolkits.axes_grid1 import AxesGrid

indices = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3), (2, 0), (2, 1), (2, 2), (2, 3), (3, 0), (3, 1), (3, 2), (3, 3)]

for mu in np.arange(0., 8.01, 1.):
	fig = plt.figure()

	grid = AxesGrid(fig, 111,
	                nrows_ncols=(4, 4),
	                axes_pad=0.025,
	                share_all=True,
	                label_mode="L",
	                cbar_location="right",
	                cbar_mode="single",
	                )
	for i, j in indices:
		a = np.array([[0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0]])
		a[i, j] = mu 
		ans = dct(dct(a, axis=0, norm='ortho'), axis=1, norm='ortho')

		pi = lambda s, a, arr : math.exp(arr[a, s])/ (math.exp(arr[0, s]) + math.exp(arr[1, s]) + math.exp(arr[2, s]) + math.exp(arr[3, s])) 

		policy = np.array([[pi(0, 0, ans), pi(0, 1, ans), pi(0, 2, ans), pi(0, 3, ans)],
					  		[pi(1, 0, ans), pi(1, 1, ans), pi(1, 2, ans), pi(1, 3, ans)],
					  		[pi(2, 0, ans), pi(2, 1, ans), pi(2, 2, ans), pi(2, 3, ans)],
					  		[pi(3, 0, ans), pi(3, 1, ans), pi(3, 2, ans), pi(3, 3, ans)]])
		im = grid[i*4 + j].imshow(policy, vmin=0, vmax=1, cmap='gray', extent = [0, 3, 0, 3])

	grid.cbar_axes[0].colorbar(im)

	for cax in grid.cbar_axes:
	    cax.toggle_label(False)

	plt.suptitle('mu = ' + str(mu))
	plt.show()
	# plt.savefig('./dct_exp_img/mu=' + str(mu) + '.png')

# # Make .gif of the images created

# import os
# import imageio

# png_dir = "./dct_exp_img/"
# images = []
# for subdir, dirs, files in os.walk(png_dir):
#     for file in files:
#         file_path = os.path.join(subdir, file)
#         if file_path.endswith(".png"):
#             images.append(imageio.imread(file_path))
# imageio.mimsave('./gif/fourByFour.gif', images, duration=2)