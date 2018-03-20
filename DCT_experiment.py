import numpy as np
import math
from scipy.fftpack import idct, dct

np.set_printoptions(precision=3, suppress=True)

print("=======Two by two (changing diagonal) ========")

for mu in np.arange(0., 8.01, 1.):
	a = np.array([[mu, 0.0],
				  [0.0, mu]]) #2.0*mu]])
	ans = dct(dct(a, axis=0, norm='ortho'), axis=1, norm='ortho')

	pi = lambda s, a, arr : math.exp(arr[a, s])/ (math.exp(arr[0, s]) + math.exp(arr[1, s])) 
	policy = np.array([[pi(0, 0, ans), pi(0, 1, ans)],
				  		[pi(1, 0, ans), pi(1, 1, ans)]]) #2.0*mu]])
	print(policy)


print("=======Three by three (changing diagonal) ========")

for mu in np.arange(0., 8.01, 1.):
	a = np.array([[mu, 0.0, 0.0],
				  [0.0, mu, 0.0],
				  [0.0, 0.0, mu]]) #2.0*mu]])
	ans = dct(dct(a, axis=0, norm='ortho'), axis=1, norm='ortho')

	pi = lambda s, a, arr : math.exp(arr[a, s])/ (math.exp(arr[0, s]) + math.exp(arr[1, s]) + math.exp(arr[2, s])) 

	policy = np.array([[pi(0, 0, ans), pi(0, 1, ans), pi(0, 2, ans)],
				  		[pi(1, 0, ans), pi(1, 1, ans), pi(1, 2, ans)],
				  		[pi(2, 0, ans), pi(2, 1, ans), pi(2, 2, ans)]]) #2.0*mu]])
	print(policy)

print("=======Three by three (changing (1, 1, 1)) ========")

for mu in np.arange(0., 8.01, 1.):
	a = np.array([[0.0, 0.0, 0.0],
				  [0.0, 0.0, 0.0],
				  [0.0, 0.0, mu]]) #2.0*mu]])
	ans = dct(dct(a, axis=0, norm='ortho'), axis=1, norm='ortho')

	pi = lambda s, a, arr : math.exp(arr[a, s])/ (math.exp(arr[0, s]) + math.exp(arr[1, s]) + math.exp(arr[2, s])) 

	policy = np.array([[pi(0, 0, ans), pi(0, 1, ans), pi(0, 2, ans)],
				  		[pi(1, 0, ans), pi(1, 1, ans), pi(1, 2, ans)],
				  		[pi(2, 0, ans), pi(2, 1, ans), pi(2, 2, ans)]]) #2.0*mu]])
	print(policy)

print("=======Four by four (changing diagonal) ========")

for mu in np.arange(0., 8.01, 1.):
	a = np.array([[mu, 0.0, 0.0, 0.0],
				  [0.0, mu, 0.0, 0.0],
				  [0.0, 0.0, mu, 0.0],
				  [0.0, 0.0, 0.0, mu]]) #2.0*mu]])
	ans = dct(dct(a, axis=0, norm='ortho'), axis=1, norm='ortho')

	pi = lambda s, a, arr : math.exp(arr[a, s])/ (math.exp(arr[0, s]) + math.exp(arr[1, s]) + math.exp(arr[2, s]) + math.exp(arr[3, s])) 

	policy = np.array([[pi(0, 0, ans), pi(0, 1, ans), pi(0, 2, ans), pi(0, 3, ans)],
				  		[pi(1, 0, ans), pi(1, 1, ans), pi(1, 2, ans), pi(1, 3, ans)],
				  		[pi(2, 0, ans), pi(2, 1, ans), pi(2, 2, ans), pi(2, 3, ans)],
				  		[pi(3, 0, ans), pi(3, 1, ans), pi(3, 2, ans), pi(3, 3, ans)]]) #2.0*mu]])
	print(policy)

print("=======Four by four (changing (1, 1, 1, 1)) ========")

for mu in np.arange(0., 8.01, 1.):
	a = np.array([[0.0, 0.0, 0.0, 0.0],
				  [0.0, 0.0, 0.0, 0.0],
				  [0.0, 0.0, 0.0, 0.0],
				  [0.0, 0.0, 0.0, mu]]) #2.0*mu]])
	ans = dct(dct(a, axis=0, norm='ortho'), axis=1, norm='ortho')

	pi = lambda s, a, arr : math.exp(arr[a, s])/ (math.exp(arr[0, s]) + math.exp(arr[1, s]) + math.exp(arr[2, s]) + math.exp(arr[3, s])) 

	policy = np.array([[pi(0, 0, ans), pi(0, 1, ans), pi(0, 2, ans), pi(0, 3, ans)],
				  		[pi(1, 0, ans), pi(1, 1, ans), pi(1, 2, ans), pi(1, 3, ans)],
				  		[pi(2, 0, ans), pi(2, 1, ans), pi(2, 2, ans), pi(2, 3, ans)],
				  		[pi(3, 0, ans), pi(3, 1, ans), pi(3, 2, ans), pi(3, 3, ans)]]) #2.0*mu]])
	print(policy)

print("=======Two by two (changing (1, 0) )========")

for mu in np.arange(0., 8.01, 1.):
	a = np.array([[0.0, 0.0],
				  [mu, 0.0]]) #2.0*mu]])
	ans = dct(dct(a, axis=0, norm='ortho'), axis=1, norm='ortho')

	pi = lambda s, a, arr : math.exp(arr[a, s])/ (math.exp(arr[0, s]) + math.exp(arr[1, s])) 
	policy = np.array([[pi(0, 0, ans), pi(0, 1, ans)],
				  		[pi(1, 0, ans), pi(1, 1, ans)]]) #2.0*mu]])
	print(policy)

print("=======Two by two (changing (0, 1) )========")

for mu in np.arange(0., 8.01, 1.):
	a = np.array([[0.0, mu],
				  [0.0, 0.0]]) #2.0*mu]])
	ans = dct(dct(a, axis=0, norm='ortho'), axis=1, norm='ortho')

	pi = lambda s, a, arr : math.exp(arr[a, s])/ (math.exp(arr[0, s]) + math.exp(arr[1, s])) 
	policy = np.array([[pi(0, 0, ans), pi(0, 1, ans)],
				  		[pi(1, 0, ans), pi(1, 1, ans)]]) #2.0*mu]])
	print(policy)

print("=======Two by two (changing (1, 0) and (1, 1))========")

for mu in np.arange(0., 8.01, 1.):
	a = np.array([[0.0, 0.0],
				  [mu, mu]]) #2.0*mu]])
	ans = dct(dct(a, axis=0, norm='ortho'), axis=1, norm='ortho')

	pi = lambda s, a, arr : math.exp(arr[a, s])/ (math.exp(arr[0, s]) + math.exp(arr[1, s])) 
	policy = np.array([[pi(0, 0, ans), pi(0, 1, ans)],
				  		[pi(1, 0, ans), pi(1, 1, ans)]]) #2.0*mu]])
	print(policy)

print("=======Two by two (changing (0, 1) and (1, 1))========")

for mu in np.arange(0., 8.01, 1.):
	a = np.array([[0.0, mu],
				  [0.0, mu]]) #2.0*mu]])
	ans = dct(dct(a, axis=0, norm='ortho'), axis=1, norm='ortho')

	pi = lambda s, a, arr : math.exp(arr[a, s])/ (math.exp(arr[0, s]) + math.exp(arr[1, s])) 
	policy = np.array([[pi(0, 0, ans), pi(0, 1, ans)],
				  		[pi(1, 0, ans), pi(1, 1, ans)]]) #2.0*mu]])
	print(policy)

print("=======Two by two (changing (1, 1))========")

for mu in np.arange(0., 8.01, 1.):
	a = np.array([[0.0, 0.0],
				  [0.0, mu]]) #2.0*mu]])
	ans = dct(dct(a, axis=0, norm='ortho'), axis=1, norm='ortho')

	pi = lambda s, a, arr : math.exp(arr[a, s])/ (math.exp(arr[0, s]) + math.exp(arr[1, s])) 
	policy = np.array([[pi(0, 0, ans), pi(0, 1, ans)],
				  		[pi(1, 0, ans), pi(1, 1, ans)]]) #2.0*mu]])
	print(policy)

print("=======Two by two (no dct)========")

for mu in np.arange(0., 8.01, 1.):
	ans = np.array([[1.0*mu, -1.0*mu],
				  [-1.0*mu, 1.0*mu]]) #2.0*mu]])

	pi = lambda s, a, arr : math.exp(arr[a, s])/ (math.exp(arr[0, s]) + math.exp(arr[1, s])) 
	policy = np.array([[pi(0, 0, ans), pi(0, 1, ans)],
				  		[pi(1, 0, ans), pi(1, 1, ans)]]) #2.0*mu]])
	print(policy)

# a = np.array([ [0.333333,  0.333333, 0.333333], \
# 			  [0.333333,  0.333333, 0.333333],
# 			  [0.333333, 0.333333, 0.333333]])
# ans = idct(idct(a, axis=0, norm='ortho'), axis=1, norm='ortho')
# print(ans)