import numpy as np

# atomic number Cu Al
Z = np.array([29., 13.])

m = .511 # GeV
E = 2.9 # GeV (elsa max energy)
p = np.sqrt(E**2 - m**2) # GeV 

beta_c = 1. # v = 1

# radiation length Cu Al
x0_rho = np.array([12.86, 24.01]) # g/cm²
x0 = np.array([1.436, 8.897]) # cm

# thickness of material
x = np.array([.07, .7]) # cm

theta = np.ndarray(shape=(2,))
first_sigma = np.ndarray(shape=(2,))
second_sigma = np.ndarray(shape=(2,))
third_sigma = np.ndarray(shape=(2,))

material = ["Copper", "Aluminium"]

# formula for multiple scattering
# first Cu second Al
for i in range(2):
    theta[i] = .0306 * 1./beta_c * 1./p * Z[i] * np.sqrt(x[i] / x0[i]) * (1. + 0.038 * np.log(x[i] / x0[i]))
    print(material[i], x[i], "[cm]")
    print("angle: ", theta[i]*180./np.pi, "[degree]")
    # 10 cm², roughly 4x4 cm² for error
    first_sigma[i] = 4./np.tan(1. * theta[i])
    second_sigma[i] = 4./np.tan(2. * theta[i])
    third_sigma[i] = 4./np.tan(3. * theta[i])
    print("first sigma: ", first_sigma[i], "[cm]")
    print("second sigma: ", second_sigma[i], "[cm]")
    print("third sigma: ", third_sigma[i], "[cm]")
    print("----------")

# questions:
# what's the max distance between the material and the detector which can be achieved at the setup?
