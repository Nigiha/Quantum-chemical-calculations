import numpy as np
import json



symbol_to_Z = {
    # 第1周期
    "H": 1, "He": 2,
    # 第2周期
    "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    # 第3周期
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    # 第4周期
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, 
    "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, 
    "Br": 35, "Kr": 36,
    # 第5周期
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, 
    "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, 
    "I": 53, "Xe": 54
}


def load_xyz(xyz_file, to_bohr=True):

    load_xyz_data=[]
    num_atoms:int
    num_eles=0

    with open(xyz_file, "r", encoding="utf-8") as f:
        lines=f.readlines()
    
    num_atoms=int(lines[0].strip())

    
    molcule_name=lines[1].strip()

    print("対象の分子 : "+molcule_name)

    electric_charge=molcule_name.count("+")-molcule_name.count("-")
    num_eles-=electric_charge

    for line in lines[2:]:
        if not line.strip():
            continue

        parts=line.split()
        symbol=parts[0]
        Z=symbol_to_Z[symbol]
        x, y, z=float(parts[1]), float(parts[2]), float(parts[3])

        num_eles+=Z

        angstrom_to_bohr=1.8897259886
        if to_bohr:
            x*=angstrom_to_bohr
            y*=angstrom_to_bohr
            z*=angstrom_to_bohr

        load_xyz_data.append({
            "symbol": symbol,
            "Z": Z,
            "coord":np.array([x, y, z])
        })

    return load_xyz_data, num_atoms, electric_charge




#二重階乗(規格化factor部分で使う)
fac2_=[0]
f_odd=1
f_even=1
for i in range(1, 5):
    if i%2==0:
        f_even*=i
        fac2_.append(f_even)
    else:
        f_odd*=i
        fac2_.append(f_odd)

def fac2(n):
    if n>0:
        return fac2_[n]
    else:
        return 1

def N_factor(a:float, l:tuple):
    l_x, l_y, l_z=l
    L=l_x+l_y+l_z

    numerartor=pow(2*a/np.pi, 0.75) * pow(4*a, L/2)
    denominator=np.sqrt(fac2(2*l_x-1)*fac2(2*l_y-1)*fac2(2*l_z-1))

    return numerartor / denominator


def load_atom_basis(json_file, atom_symbol:str):
    
    extracted_shells=[]

    with open(json_file, "r", encoding="utf-8") as f:
        basis_data=json.load(f)

    Z_str=str(symbol_to_Z[atom_symbol])

    if Z_str not in basis_data["elements"]:
        raise ValueError(f"{atom_symbol} の基底関数データが見つかりません。")
    
    shells=basis_data["elements"][Z_str]["electron_shells"]

    for shell in shells:
        
        angular_momentum=shell["angular_momentum"]
        
        for i in range(len(angular_momentum)):
            
            l=angular_momentum[i]
            for l_x in range(0, l+1):
                for l_y in range(0, l-l_x+1):
                    l_z=l-l_x-l_y
                    
                    angular_momentum_vector=(l_x, l_y, l_z)
                    exponents=[float(e) for e in shell["exponents"]]
                    coefficients=[float(c) for c in shell["coefficients"][i]]

                    normalization_factors=[N_factor(float(e), angular_momentum_vector) for e in shell["exponents"]]

                    extracted_shells.append({
                        "angular_momentum": angular_momentum_vector,
                        "exponents": exponents,
                        "coefficients": coefficients,
                        "normalization_factors": normalization_factors,
                    })

    return extracted_shells

# print(load_atom_basis("STO-3G.json", "O"))
