import numpy as np

# atomic number Cu Al
Z = np.array([29., 13.])

e = 1.602e-19 # C
c = 2.998e8 # m/s
m = .000511 # GeV
m_si = 9.109e-31 # kg
E = 2.9 # GeV (elsa max energy)
E_si = (E * 1e-9) * e
#v = np.sqrt(-1 * (1 + E_si / (m_si * c**2)) + 1) * c
#print("vel: ",v)
p = np.sqrt(E**2 - m**2) # GeV 

beta_c = 1. # v = 1

# radiation length Cu Al
x0_rho = np.array([12.86, 24.01]) # g/cm²
x0 = np.array([1.436, 8.897]) # cm

# thickness of material
x = np.array([10, 7000]) # cm

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
    print("angle: ", theta[i], "[degree]")
    # 10 cm², roughly 4x4 cm² for error
    first_sigma[i] = 4./np.tan(1. * theta[i]*np.pi/180)
    second_sigma[i] = 4./np.tan(2. * theta[i]*np.pi/180)
    third_sigma[i] = 4./np.tan(3. * theta[i]*np.pi/180)
    print("first sigma: ", first_sigma[i], "[cm]")
    print("second sigma: ", second_sigma[i], "[cm]")
    print("third sigma: ", third_sigma[i], "[cm]")
    print("----------")

''' testing
# cummulated
theta_c = np.ndarray(shape=(3,2))
x_c = np.array([[.01, .1], [.02, .2], [.04, .4]])
third_sigma_c = np.ndarray(shape=(2,3))
resolution = np.ndarray(shape=(2,))

test = np.array([[0,1],[2,3],[4,5]])
print("asdf", test[0][1])

for j in range(3):
    for i in range(2):
        theta_c[j][i] = .0306 * 1./beta_c * 1./p * Z[i] * np.sqrt(x_c[j][i] / x0[i]) * (1. + 0.038 * np.log(x_c[j][i] / x0[i]))

theta_med = np.sqrt(theta_c[i][j].sum())
print(theta_med)

resolution[0] = third_sigma_[0][0]*(np.tan(theta_c[0][0]))


# questions:
# what's the max distance between the material and the detector which can be achieved at the setup?
'''
