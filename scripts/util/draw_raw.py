#f:\\temp\\Gas_smoothlength_100_100_100_float.vdb
#f:\\temp\\Gas_BFieldz_100_100_100_float.vdb
#f:\\temp\\Gas_BFieldz_100_100_100_float.raw

import numpy as np
import matplotlib.pyplot as plt

def read_vdb(filename, dim):
    import pyopenvdb as vdb

    grid = vdb.read(filename, gridname='density')
    array = np.ndarray((dim, dim, dim), np.float32)
    grid.copyToArray(array, ijk=(0, 0, 0))

    array = np.transpose(array, axes=(2, 1, 0))  # New shape: (z, y, x)

    return array

def read_raw(filename, dim):
    array = np.fromfile(filename, dtype=np.float32)
    array = array.reshape((dim, dim, dim))    

    return array

def print_info2(array1, array2):
    print(array1.shape, array2.shape)
    print(np.sqrt(np.sum((array1 - array2) ** 2)))

    # Compute absolute differences
    #differences = np.abs(array1 - array2).flatten()
    # top_n = 1000
    # largest_differences = np.sort(differences)[-top_n:][::-1]
    # print(largest_differences)

    print("array1", np.min(array1), np.max(array1))
    print("array2", np.min(array2), np.max(array2))

def print_info1(array1):
    print(array1.shape)
    print("array1", np.min(array1), np.max(array1))

def draw_slice_ramp(array1, s):
    #vmin = np.min(array1) / 10000
    #vmax = np.max(array1) / 10000
    
    # Extract the slice at index 50 along the first axis
    slice_50 = array1[s, :, :]

    # Plot the slice using a color ramp
    plt.figure(figsize=(10, 10))
    plt.imshow(slice_50, cmap='viridis') #, vmin=vmin, vmax=vmax)  # You can use other colormaps like 'plasma', 'inferno', etc.
    plt.colorbar(label="Intensity")  # Add a colorbar to show the intensity scale
    plt.title(f"Slice at Index {s} (x={s}) with Color Ramp")
    plt.xlabel("Y-axis")
    plt.ylabel("Z-axis")
    plt.show()

def draw_radius(array1, max_radius):
    # Extract the slice at index 50 along the first axis
    slice_50 = array1[50, :, :]

    # Create a figure
    fig, ax = plt.subplots(figsize=(10, 10))

    # Iterate through the 2D slice and plot circles for non-zero values
    for y in range(slice_50.shape[0]):  # Iterate over rows
        for x in range(slice_50.shape[1]):  # Iterate over columns
            value = slice_50[y, x]
            if value > 0:  # Draw a circle only for non-zero values
                circle = plt.Circle((x, y), radius=value/max_radius, color='blue', alpha=0.5)
                ax.add_artist(circle)

    # Set limits to fit the grid
    ax.set_xlim(0, slice_50.shape[1])
    ax.set_ylim(0, slice_50.shape[0])
    ax.set_aspect('equal')  # Keep the aspect ratio square

    # Invert the y-axis to match the matrix layout
    ax.invert_yaxis()

    # Add labels and title
    plt.title("Circles Representing Non-Zero Values (Radius = Value)")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")

    plt.show()

def voxelization2D():
    # Define parameters
    radius = 1.0  # Radius of the sphere
    origin_value = 1.0  # Original value of the sphere
    grid_size = 9  # Number of pixels along one dimension (3x3 grid)
    sigma = 0.5  # Standard deviation for the Gaussian kernel

    # Create a grid of coordinates
    x = np.linspace(-1.5 * radius, 1.5 * radius, grid_size)
    y = np.linspace(-1.5 * radius, 1.5 * radius, grid_size)
    xx, yy = np.meshgrid(x, y)

    # Calculate distances from the center (0, 0)
    distances = np.sqrt(xx**2 + yy**2)

    # Apply Gaussian kernel
    gaussian_kernel = np.exp(-distances**2 / (2 * sigma**2))

    # Normalize the kernel so the sum matches the original sphere value
    voxel_values = gaussian_kernel / np.sum(gaussian_kernel) * origin_value

    # Display the voxelized grid
    plt.imshow(voxel_values, cmap='hot', interpolation='nearest', origin='lower')
    plt.colorbar(label='Pixel Value')
    plt.title('Voxelized 2D Sphere with Gaussian Kernel')
    plt.xticks(range(grid_size), labels=np.round(x, 2))
    plt.yticks(range(grid_size), labels=np.round(y, 2))
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.show()

    # Print the voxel values and their sum
    #print("Voxel values:\n", voxel_values)
    print("Sum of voxel values:", np.sum(voxel_values))

def gaussian_kernel_3d(size: int, sigma: float):
    """
    Generate a 3D Gaussian kernel.
    
    Parameters:
    - size: int, the size of the kernel along one dimension (must be odd).
    - sigma: float, the standard deviation of the Gaussian distribution.
    
    Returns:
    - kernel: 3D numpy array of shape (size, size, size).
    """
    ax = np.linspace(-(size // 2), size // 2, size)
    x, y, z = np.meshgrid(ax, ax, ax, indexing='ij')
    kernel = np.exp(-(x**2 + y**2 + z**2) / (2 * sigma**2))
    kernel /= np.sum(kernel)
    return kernel

def gauss3d():
    # Create a 3D Gaussian kernel with size 3x3x3 and sigma 1
    kernel_size = 9
    sigma = 1.0
    gaussian_kernel = gaussian_kernel_3d(kernel_size, sigma)

    # Verify that the sum of all elements equals 1
    sum_of_elements = np.sum(gaussian_kernel)

    print("3D Gaussian Kernel:")
    #print(gaussian_kernel)
    print("\nSum of all elements in the kernel:", sum_of_elements)

    # Check if the sum is approximately 1
    if np.isclose(sum_of_elements, 1.0):
        print("\nThe sum of all elements in the Gaussian kernel is approximately 1.")
    else:
        print("\nThe sum of all elements in the Gaussian kernel is not approximately 1.")

def test_gaus3d():
    import numpy as np
    import math

    # Parameters
    px, py, pz = 4.0, 4.0, 4.0
    iradiusx, iradiusy, iradiusz = 2, 2, 2

    # Initialize the kernels array and sum
    kernels_size = (iradiusx + 1) * (iradiusy + 1) * (iradiusz + 1) * 8
    kernels = np.zeros(kernels_size)
    kernels_sum = 0.0

    # Populate the kernels array
    for sx in range(int(px - iradiusx), int(px + iradiusx) + 1):
        for sy in range(int(py - iradiusy), int(py + iradiusy) + 1):
            for sz in range(int(pz - iradiusz), int(pz + iradiusz) + 1):
                dx = sx - px
                dy = sy - py
                dz = sz - pz

                sigma = 1.0
                W = math.exp(-1.0 * (dx**2 + dy**2 + dz**2) / (2 * sigma**2))  # Gaussian kernel

                density = W
                norm = W

                index = (
                    iradiusx * iradiusy * (sz - pz + iradiusz) +
                    iradiusx * (sy - py + iradiusy) +
                    (sx - px + iradiusx)
                )
                kernels[int(index)] = density
                kernels_sum += norm

    # Check the sum
    density_check = 0.0
    for sx in range(int(px - iradiusx), int(px + iradiusx) + 1):
        for sy in range(int(py - iradiusy), int(py + iradiusy) + 1):
            for sz in range(int(pz - iradiusz), int(pz + iradiusz) + 1):
                index = (
                    iradiusx * iradiusy * (sz - pz + iradiusz) +
                    iradiusx * (sy - py + iradiusy) +
                    (sx - px + iradiusx)
                )
                W = kernels[int(index)]
                density_check += W / kernels_sum

    # Verify the result
    if abs(density_check - 1.0) > 1e-5:
        print(f"Error: check sum is not one: {density_check}")
    else:
        print("Density check passed. Sum is approximately 1.0")

def test2():
    import numpy as np
    import math

    # Parameters
    px, py, pz = 4.0, 4.0, 4.0
    iradiusx, iradiusy, iradiusz = 2, 2, 2

    # Determine the size of the kernel based on the radius
    kernel_dim_x = 2 * iradiusx + 1
    kernel_dim_y = 2 * iradiusy + 1
    kernel_dim_z = 2 * iradiusz + 1

    # Initialize the kernels array
    kernels = np.zeros((kernel_dim_x, kernel_dim_y, kernel_dim_z))
    kernels_sum = 0.0

    # Populate the kernels array
    for sx in range(kernel_dim_x):
        for sy in range(kernel_dim_y):
            for sz in range(kernel_dim_z):
                dx = sx - iradiusx
                dy = sy - iradiusy
                dz = sz - iradiusz

                sigma = 1.0
                W = math.exp(-1.0 * (dx**2 + dy**2 + dz**2) / (2 * sigma**2))  # Gaussian kernel

                kernels[sx, sy, sz] = W
                kernels_sum += W

    # Normalize the kernel
    kernels /= kernels_sum

    # Check the sum
    density_check = np.sum(kernels)

    # Verify the result
    if abs(density_check - 1.0) > 1e-5:
        print(f"Error: check sum is not one: {density_check}")
    else:
        print("Density check passed. Sum is approximately 1.0")

def test3():
    import numpy as np
    import math

    # Parameters
    px, py, pz = 4.0, 4.0, 4.0  # Center position
    iradiusx, iradiusy, iradiusz = 2, 2, 2  # Radii

    # Kernel dimensions
    kernel_dim_x = 2 * iradiusx + 1
    kernel_dim_y = 2 * iradiusy + 1
    kernel_dim_z = 2 * iradiusz + 1

    # Initialize the kernel and sum
    kernels = np.zeros((kernel_dim_x, kernel_dim_y, kernel_dim_z))
    kernels_sum = 0.0

    # Populate the kernel array
    for sx in range(kernel_dim_x):
        for sy in range(kernel_dim_y):
            for sz in range(kernel_dim_z):
                dx = sx - iradiusx
                dy = sy - iradiusy
                dz = sz - iradiusz

                sigma = 1.0
                # Gaussian kernel formula
                W = math.exp(-1.0 * (dx**2 + dy**2 + dz**2) / (2 * sigma**2))

                kernels[sx, sy, sz] = W
                kernels_sum += W

    # Normalize the kernel so the sum equals 1
    kernels /= kernels_sum

    # Verify the normalization
    density_check = np.sum(kernels)

    # Print results
    if abs(density_check - 1.0) > 1e-6:
        print(f"Error: check sum is not one: {density_check}")
    else:
        print("Density check passed. Sum is approximately 1.0")

    # Optional: Print the kernel values
    #print("3D Gaussian Kernel:")
    #print(kernels)


array1 = read_raw("f:\\temp\\Gas_B_50_50_50_float.raw", 50)
#array2 = read_vdb("f:\\temp\\Gas_B.vdb", 50)
#print_info2(array1, array2)
#for s in range(0, 50, 10):
#    draw_slice_ramp(array1,s)
draw_slice_ramp(array1,25)

# array3 = read_vdb("f:\\temp\\Gas_smoothlength_100_100_100_float.vdb", 100)
# print_info1(array3)
# print_info1(array3[50, :, :])
#draw_radius(array3, 1490.9967)

#voxelization2D()
#gauss3d()
#test_gaus3d()
#test2()
#test3()