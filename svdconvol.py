import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import svd


# load png image as numpy matrix at given user path
def load_image(path):
    from PIL import Image
    img = Image.open(path)
    img.load()
    data = np.asarray(img, dtype="float32")
    return data
    
# smoothstep function
def smoothstep(min, max, x):
    if x < min:
        return 0.0
    if x > max:
        return 1.0
    x = (x - min) / (max - min)
    return x * x * (3 - 2 * x)


matrix = load_image("C:/Users/micha/Downloads/Capture.png")
matrix = matrix[:,:,0]

#remove matrix rows and columns until its square
while matrix.shape[0] != matrix.shape[1]:
    if matrix.shape[0] > matrix.shape[1]:
        matrix = np.delete(matrix, -1, 0)
    else:
        matrix = np.delete(matrix, -1, 1)

N = matrix.shape[0]
Nc = N//2

# define the N*N matrix with a Gaussian in it
#N = 91
#Nc = 45
#matrix = np.zeros((N, N))
#for i in range(N):
#    for j in range(N):
#        R = np.sqrt((i - Nc) ** 2 + (j - Nc) ** 2)
#        # gaussian line in x direction
#        xspike = 0.5 * np.exp(-((i - Nc) / 10.0) ** 2)
#        # gaussian line in y direction
#        yspike = 0.5 * np.exp(-((j - Nc) / 10.0) ** 2)
#        kernel = xspike + yspike + (R**2 + 8.0)**(-1) 
#        matrix[i, j] = (1.0 - smoothstep(N*0.45, N * 0.5, R))

# normalize the matrix, so that the sum of all elements is 1
matrix /= np.sum(matrix)

# calculate the SVD decomposition
U, s, V = svd(matrix)

num_ranks = 2

# print out the U matrix columns in a form of a glsl array
print("//The first {} columns of the U matrix:".format(num_ranks))

for i in range(num_ranks):
    print("float U{}[{}] = float[](".format(i + 1, N))
    array = ""
    for j in range(N-1):
        array += "    {}, ".format(U[j, i])
    array += "    {} ".format(U[N-1, i])
    print(array)
    print(");")

# print out the V matrix rows in a form of a glsl array

print("//The first {} rows of the V matrix:".format(num_ranks))
for i in range(num_ranks):
    print("float V{}[{}] = float[](".format(i + 1, N))
    array = ""
    for j in range(N-1):
        array += "    {}, ".format(V[i, j])
    array += "    {} ".format(V[i, N-1])
    print(array)
    print(");")

# print out the singular values in a form of a glsl array
print("//The first {} singular values:".format(num_ranks))
print("float S[{}] = float[](".format(num_ranks))
array = ""
for i in range(num_ranks-1):
    array += "    {}, ".format(s[i])
array += "    {} ".format(s[num_ranks-1])
print(array)
print(");")

#print out Nc in glsl form
print("//center of the convolution")
print("int Nc = {};".format(Nc))


# reconstruct the matrix using 2 first ranks
matrix_reconstructed = np.zeros((N, N))
for i in range(num_ranks):
    matrix_reconstructed += s[i] * np.outer(U[:, i], V[i, :])

# plot the original and reconstructed matrices
fig, axs = plt.subplots(1, 2)
axs[0].imshow(matrix, cmap='gray')
axs[0].set_title("Original matrix")
axs[1].imshow(matrix_reconstructed, cmap='gray')
axs[1].set_title("Reconstructed matrix (using {} ranks)".format(num_ranks))
plt.show()

