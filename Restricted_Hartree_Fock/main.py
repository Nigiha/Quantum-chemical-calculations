import numpy as np
from scipy.integrate import quad

#==========step1:対象分子の設定==========


target_molecule_file="HeH.xyz"
base_function_file="STO-3G.json"


# 元素記号から原子番号への変換辞書リスト
symbol_to_Z = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10
}



#核電荷{Z}と核座標{R=(x,y,z)}の受け取り

angstrom_to_bohr=1.8897259886


def load_xyz(xyz_file, to_bohr=True):
    molecule_data=[]

    with open(xyz_file, "r", encoding="utf-8") as f:
        lines=f.readlines()
    
    num_atoms=int(lines[0].strip())

    for line in lines[2:]:
        if not line.strip():
            continue

        parts=line.split()
        symbol=parts[0]
        x, y, z=float(parts[1]), float(parts[2]), float(parts[3])

        if to_bohr:
            x*=angstrom_to_bohr
            y*=angstrom_to_bohr
            z*=angstrom_to_bohr

        molecule_data.append({
            "symbol": symbol,
            "coord":np.array([x, y, z])
        })

    return molecule_data, num_atoms
    
molecule, N=load_xyz(target_molecule_file, to_bohr=True)



#基底関数の受け取り(STO-3G.jsonからGTOの指数と係数のデータを抜粋)
import json
def load_atom_basis(json_file, atom):
    with open(json_file, "r", encoding="utf-8") as f:
        basis_data=json.load(f)

    Z_str=str(symbol_to_Z[atom])

    if Z_str not in basis_data["elements"]:
        raise ValueError(f"{atom} の基底関数データが見つかりません。")
    
    shells=basis_data["elements"][Z_str]["electron_shells"]

    extracted_shells=[]
    for shell in shells:
        exponents=[float(e) for e in shell["exponents"]]
        coefficients=[float(c) for c in shell["coefficients"][0]] #今回はs軌道のみを考える

        extracted_shells.append({
            "exponents": exponents,
            "coefficients":coefficients
        })

    return extracted_shells



#分子全体の基底関数リストを作る
def build_molecule_basis(molecule, json_file):
    molecule_basis=[]

    for atom in molecule:
        symbol=atom["symbol"]
        xyz=atom["coord"]
        Z=symbol_to_Z[symbol]

        atom_shells = load_atom_basis(json_file, symbol)

        for shell in atom_shells:
            molecule_basis.append({
                "atom_symbol": symbol,
                "Z": Z,
                "center": xyz,
                "exponents": shell["exponents"],
                "coefficients": shell["coefficients"]
            })
    
    return molecule_basis

molecule_basis=build_molecule_basis(molecule, base_function_file)
# print(molecule_basis)



#核間反発エネルギーV_nucの計算
# def V_nuc_matrix(molecule):

#     V=np.zeros((N, N))
    
#     for i in range(N):
#         for j in range(i+1, N):
#             r_i, r_j=np.array(molecule[i]["coord"]), np.array(molecule[j]["coord"])
#             Z_i, Z_j=symbol_to_Z[molecule[i]["symbol"]], symbol_to_Z[molecule[j]["symbol"]]
#             V[i][j]=Z_i*Z_j/(np.linalg.norm(r_i-r_j))
#             V[j][i]=Z_i*Z_j/(np.linalg.norm(r_i-r_j))
    
#     return V

def V_nuc_schalar(molecule):

    V=0
    
    for i in range(N):
        for j in range(i+1, N):
            r_i, r_j=np.array(molecule[i]["coord"]), np.array(molecule[j]["coord"])
            Z_i, Z_j=symbol_to_Z[molecule[i]["symbol"]], symbol_to_Z[molecule[j]["symbol"]]
            V+=Z_i*Z_j/(np.linalg.norm(r_i-r_j))
    
    return V


K=len(molecule_basis) #基底関数の数

def Nm_c(a): #GTOの規格化定数
    return (2*a/np.pi)**(3/4)

#==========step2:正準直交化の実行==========


def S_prim(a:float, b:float, A:np.ndarray, B:np.ndarray): #原始GTO同士の重なり積分　←　多中心積分の公式とガウス積分から
        return (np.pi/(a+b))**(3/2) * np.exp(-a*b/(a+b)*np.linalg.norm(A-B)**2)


# 重なり行列Sの計算(GTOをs型に限れば、解析計算は容易)
def S(molecule_basis):
    S=np.zeros((K, K))

    for m in range(K):
        for n in range(m, K):
            alpha_m, alpha_n=molecule_basis[m]["exponents"], molecule_basis[n]["exponents"]
            c_m, c_n=molecule_basis[m]["coefficients"], molecule_basis[n]["coefficients"]

            A_m=molecule_basis[m]["center"]
            A_n=molecule_basis[n]["center"]

            s=0

            for i in range(len(alpha_m)):
                for j in range(len(alpha_n)):
                    s+=c_m[i]*c_n[j]*Nm_c(alpha_m[i])*Nm_c(alpha_n[j])*S_prim(alpha_m[i], alpha_n[j], A_m, A_n)
            
            S[n][m]=s
            S[m][n]=s
    
    return S



#Sの対角化
val_S, vec_S=np.linalg.eigh(S(molecule_basis))

val_S_inv_sqrt = 1.0 / np.sqrt(val_S)
s_inv_sqrt=np.diag(val_S_inv_sqrt, 0)
U=vec_S


#変換行列Xの計算(正準直交化)                 
X=U@s_inv_sqrt





#==========step3==========


def T_prim(a:float, b:float, A:np.ndarray, B:np.ndarray):
    return 


#運動エネルギー項Tの計算
def T(molecule_basis):
    T=np.zeros((K, K))

    



#核引力項Vの計算
