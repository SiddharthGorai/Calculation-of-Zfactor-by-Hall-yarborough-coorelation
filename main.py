# this program calcaulates z-factor of a real gas by using standing's correlation and hall-yarborough method.

import numpy as np
import csv

# Initial Well data
P = 7111709.0 # Pressure in Pa
T = 338.15 # Temp in k
mol_wt_air = 28.96 #g/mol

#importing data
comp_data_file = open(r"composition_data.csv","r")
comp_csv_reader = csv.reader(comp_data_file)

mol_wt = [] # g/mol
yi_data = [] # mol fraction
Tci_data = [] # oR
Pci_data = [] # psia

for i in comp_csv_reader:
    mol_wt.append(i[1])
    yi_data.append(i[2])
    Tci_data.append(i[3])
    Pci_data.append(i[4])

mol_wt.pop(0)
yi_data.pop(0)
Tci_data.pop(0)
Pci_data.pop(0)

mol_wt = np.array(mol_wt)
yi_data = np.array(yi_data)
Tci_data =  np.array(Tci_data)
Pci_data =  np.array(Pci_data)

mol_wt = mol_wt.astype(float)
yi_data = yi_data.astype(float)
Tci_data =  Tci_data.astype(float)
Pci_data =  Pci_data.astype(float)

SG = np.sum(mol_wt*yi_data)/mol_wt_air

def calculate_Ppr_Tpr(method):
    # Calculating Psuedo Reduced pressure and temp 
   
    if method == 1:
        # Method 1
        Ppc = np.sum(yi_data*Pci_data)
        Tpc = np.sum(yi_data*Tci_data)

    if method == 2: 

        # Method 2 Sutton's Correlation
        Ppc = 756.8 - 131*SG - 3.6*(SG**2)
        Tpc = 169.2 + 349.5*SG - 74*(SG**2)


    if method == 3:
        # Method 3 Standing's Correlation

        if SG < 0.75:
            Ppc = 667 + 15*SG - 37.5*(SG**2)
            Tpc = 168 + 325*SG - 12.5*(SG**2)
        else:
            Ppc = 706 + 51.7*SG - 11.1*(SG**2)
            Tpc = 187 + 330*SG - 71.5*(SG**2)
    
    Ppc = Ppc*100000/14.7  # conversion of Psia to Pa
    Tpc = ((Tpc-460-32)*(5/9)) + 273 # conversion of oR to K

    Ppr = P/Ppc
    Tpr = T/Tpc
    return (Ppr,Tpr)


method = 1 # Calculating Ppr and TPr

Ppr = calculate_Ppr_Tpr(method)[0]
Tpr = calculate_Ppr_Tpr(method)[1]

print(Ppr,Tpr)
# Hall-Yarborough Parameters
tr = (1/Tpr)
A = 0.06125*np.power(np.e,-1.2*(1-tr**2))
B = tr*(14.96 - 9.76*tr + 4.58*(tr**2))
C = tr*(90.7 - 242.2*tr + 42.4*(tr**2))
D = 2.18 + 2.82*tr
y = 0.01
iterations = 3
# y1 = y - fy/dfy



for i in range(iterations):
    fy = ((y + y**2 + y**3 - y**4)/((1-y)**3)) - A*Ppr - B*(y**2) + C*(y**D)
    dfy = ((1 + 4*y + 4*(y**2) - 4*(y**3) + y**4 )/((1 - y)**4))  - 2*B*y + C*D*(y**(D-1))
    y = y - fy/dfy
    z = A*Ppr/y
print(z)