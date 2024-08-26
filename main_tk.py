import numpy as np
import csv
from tkinter import *
from tkinter import ttk
from tkinter import messagebox
import os
import matplotlib.pyplot as plt

def show_csv_files():
    files = []
    for i in os.listdir():
        if i.endswith(".csv"):
            files.append(i)
    return files

def show_data(event):
    for i in csv_table_labels:
        i.destroy()
    csv_table_labels.clear()
    data = []
    df = data_combobox.get()
    csv_data_file = open(df,'r')
    csv_data_reader = csv.reader(csv_data_file)
    for row in csv_data_reader:
        data.append(row)
    if data:
    
        for row_index, row in enumerate(data):
            for col_index, value in enumerate(row):
                label = ttk.Label(win, text=value)
                label.grid(row=row_index + 7, column=col_index, padx=5, pady=5)
                csv_table_labels.append(label)



def calculate_Ppr_Tpr(method,yi_data,Pci_data,Tci_data,SG,Pressure,Temp):
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

    Ppr = Pressure/Ppc
    Tpr = Temp/Tpc
    # print(Pressure, Ppc,Ppr,Tpr)
    return (Ppr,Tpr)

def show_graph():
        
        pressure_file_name = pressure_combobox.get()
        presssure_data_file = open(pressure_file_name,"r")
        pressure_csv_reader = csv.reader(presssure_data_file)

        pressure = []
        temp = []
        z = []

        for i in pressure_csv_reader:
            pressure.append(i[0])
            temp.append(i[1])
        pressure.pop(0)
        temp.pop(0)

        pressure = np.array(pressure)
        temp = np.array(temp)

        pressure = pressure.astype(float)
        temp = temp.astype(float)

    
        # print(pressure,temp)
        for p,t in zip(pressure,temp):
            # print(i,j)
            zi = calculate_z(p,t)
            z.append(zi)
        # print(z)
        z = np.array(z)

        plt.scatter(pressure,z, marker='o')
        plt.ylabel("z-factor")
        plt.xlabel("Pressure (pa)")
        plt.show()

                



def calculate_z(P,T):
    try:
        comp_file_name = data_combobox.get()
        if(comp_file_name.endswith(".csv") ):
            comp_data_file = open(comp_file_name,"r")
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

            try:
                mol_wt_air = float(molwt_entry.get())
                SG = np.sum(mol_wt*yi_data)/mol_wt_air
            except:
                messagebox.showerror("Error","Please enter valid Specific Gravity.")

            method = method_combobox.get()
            if method == "Composition":
                method = 1
            elif method == "Sutton's Correlation":
                method = 2
            elif method == "Standing's Correlation":
                method = 3
            else:
                messagebox.showerror("Error","Select valid method.")

            # P = pressure_entry.get()
            # T = temp_entry.get()

            if(not (P == "" and T == "")):
                try:
                    P = float(P)
                    T = float(T)
                    Ppr = calculate_Ppr_Tpr(method,yi_data,Pci_data,Tci_data,SG,P,T)[0]
                    Tpr = calculate_Ppr_Tpr(method,yi_data,Pci_data,Tci_data,SG,P,T)[1]
                    # Hall-Yarborough Parameters
                    tr = (1/Tpr)
                    A = 0.06125*np.power(np.e,-1.2*((1-tr)**2))*tr
                    B = tr*(14.76 - 9.76*tr + 4.58*(tr**2))
                    C = tr*(90.7 - 242.2*tr + 42.4*(tr**2))
                    D = 2.18 + 2.82*tr
                    y = float(y_entry.get())
                    while True:
                        fy = ((y + y**2 + y**3 - y**4)/((1-y)**3)) - A*Ppr - B*(y**2) + C*(y**D)
                        dfy = ((1 + 4*y + 4*(y**2) - 4*(y**3) + y**4 )/((1 - y)**4))  - 2*B*y + C*D*(y**(D-1))
                        zi = A*Ppr/y
                        y = y - fy/dfy
                        z = A*Ppr/y
                        # print(fy)
                        if(abs(zi - z) < (10**-5)):
                            break
                    mol_wt_air_in_kgpermol = mol_wt_air / 1000
                    
                    rho = (SG*mol_wt_air_in_kgpermol*P)/(8.314*T*z) # kg/m3
                    rho = rho / 1000 # g/cc
                    rho_label.config(text=" Value of Gas density is: {} g/cm3".format(rho))
                    z_label.config(text=" Value of Z is: {}".format(z))

                    return z
                    
                except:
                    messagebox.showerror("Error", "Please enter valid data")
        
        else:
            messagebox.showerror("Error","Select valid composition data.")
    except Exception as e:
        messagebox.showerror("Error",e)


win = Tk()
win.title("Z-Factor Calculator")
win.resizable(False, False)

csv_table_labels = []

head_label = Label(win,text="Z-Factor Calculator using Hall-Yarborough Method",font=("Arial",10,"bold"))
head_label.grid(row=0,column=0,columnspan=4)

pressure_label =  Label(win,text="Pressure (Pa)")
pressure_label.grid(row=1,column=0 )

pressure_entry = Entry(win)
pressure_entry.insert(0,7111709)
pressure_entry.grid(row=2,column=0, padx=10)

molwt_label =  Label(win,text="Molecular Wt of air (g/mol)")
molwt_label.grid(row=4,column=0)

molwt_entry = Entry(win)
molwt_entry.insert(0,28.96)
molwt_entry.grid(row=5,column=0, padx=10)

temp_label =  Label(win,text="Temperature (K)")
temp_label.grid(row=1,column=1)

temp_entry = Entry(win)
temp_entry.insert(0,338.15)
temp_entry.grid(row=2,column=1 )

y_label = Label(win,text="Reduced density parameter (Yo)")
y_label.grid(row=4,column=1)

y_entry = Entry(win)
y_entry.insert(0,0.01)
y_entry.grid(row=5,column=1 )

method_label =  Label(win,text="Method for Ppr and Tpr Calculation")
method_label.grid(row=1,column=2 )
method_options = ["Composition", "Sutton's Correlation", "Standing's Correlation"]
method_combobox = ttk.Combobox(win, values=method_options, state="readonly")
method_combobox.set("Select any one")  
method_combobox.grid(row=2, column=2 )

# pressure_label =  Label(win,text="Select pressure data")
# pressure_label.grid(row=1,column=3 )
# pressure_options = show_csv_files()
# pressure_combobox = ttk.Combobox(win, values=pressure_options, state="readonly")
# pressure_combobox.set("Select any one")  
# pressure_combobox.grid(row=2, column=3 )

data_label =  Label(win,text="Select CSV datafile")
data_label.grid(row=4,column=2 )
data_options = show_csv_files()
data_combobox = ttk.Combobox(win, values=data_options, state="readonly")
data_combobox.set("Select any one") 
data_combobox.bind("<<ComboboxSelected>>",show_data)
data_combobox.grid(row=5, column=2)

empty_label = Label(win)
empty_label.grid(row=6, column = 0)

P = pressure_entry.get()
T = temp_entry.get()

start_button = Button(win,text="Calculate Z",command=lambda: calculate_z(P,T))
start_button.grid(row=5,column=3,padx=20)

rho_label =  Label(win)
rho_label.grid(row=2,column=4)

z_label =  Label(win)
z_label.grid(row=3,column=4)

# graph_button = Button(win,text="Show Z/P Graph",command=show_graph)
# graph_button.grid(row=5,column=4,padx=20)

win.mainloop()
