import numpy as np

# atomic number Cu Al
Z = np.array([29., 13.])

m = .000511 # GeV
E = 2.9 # GeV (elsa max energy)
p = np.sqrt(E**2 - m**2) # GeV 

beta_c = 1. # v = 1

# radiation length Cu Al
x0_rho = np.array([12.86, 24.01]) # g/cm²
x0 = np.array([1.436, 8.897]) # cm

# thickness of material
x = np.array([2.96, 18.65]) # cm

# distance target detector
d = 40 # cm
#print("distance: ", d, "[cm]")

theta = np.ndarray(shape=(2,))
distance_target_detector = np.ndarray(shape=(2,))

material = ["Copper", "Aluminium"]

# formula for multiple scattering
# first Cu second Al
for i in range(2):
    theta[i] = .0136 * 1./beta_c * 1./p * 1 * np.sqrt(x[i] / x0[i]) * (1. + 0.038 * np.log(x[i] / x0[i]))

def multipleScattering(thickness, radLength, Z):
    return .0136 * (1./beta_c) * (1./p) * 1 * np.sqrt((thickness/radLength)) * (1.0 + 0.038 * np.log(thickness/radLength))

thicknessArray = np.array([
        5, 9.93, 21, 27.39,
        0.7, 1.48, 2.96, 4.36
        ])
materialArray = np.array([
        1,1,1,1,
        0,0,0,0
        ])

for i in range(0,8):
    theta = multipleScattering(
            thicknessArray[i], 
            x0[materialArray[i]], 
            Z[materialArray[i]]
            )
    print(material[materialArray[i]], thicknessArray[i], "[cm]")

    print("angle: ", theta, "[rad]")
    print("angle: ", theta*180/np.pi, "[degree]")
    
    print("-------------------------")
