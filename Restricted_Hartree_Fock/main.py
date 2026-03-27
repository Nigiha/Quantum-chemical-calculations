#==========step1:対象分子の設定==========

target_molecule="HeH"
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
from molecule import molecule_structure
molecule=molecule_structure[target_molecule]

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

basis=build_molecule_basis(molecule, "STO-3G.json")
print(basis)























import numpy as np
from scipy.integrate import quad

# #核間反発エネルギーV_nucの計算
# def V_nuc(N, mol_str):
#     Z, x, y, z=mol_str[0,:], mol_str[1,:], mol_str[2,:], mol_str[3,:]
#     V=np.zeros((N, N))
#     for i in range(N):
#         for j in range(i+1, N):
#             r_x=x[i]-x[j]
#             r_y=y[i]-y[j]
#             r_z=z[i]-z[j]
#             r=np.sqrt(r_x**2+r_y**2+r_z**2)
#             if r==0:
#                 print("エラー：")
#                 exit()
#             V[i][j]=Z[i]*Z[j]/r
#             V[i][j]=Z[i]*Z[j]/r
#     return V

# #==========step2:正準直交化の実行==========
# # 重なり行列Sの計算

# #Sの対角化

# #変換行列Xの計算(正準直交化)                 

# #==========step3:コアハミルトニアン行列の計算==========

# #運動エネルギー項Tの計算
# def T(N, a_i):
#     def GTO(a):
#         return lambda r:(2*a/np.pi)**(3/4)*np.exp(-a*r**2)

#     def laplace_GTO(a):
#         return lambda r:(2*a/np.pi)**(3/4) * 2*a(2*a*r**2-3)*np.exp(-a*r**2)

#     J=lambda r:4*np.pi*r**2
    
#     T=np.zeros((N, N))
#     for i in range(N):
#         for j in range(N):
#             T[i][j]=quad(GTO(a_i[i])*(1/2)*laplace_GTO(a_i[j]) * J(), -np.inf, np.inf)
#     return T

# #核引力項の計算
# def V_ne(N, a_i):
#     def xi(a):
#         return lambda x, y, z:(2*a/np.pi)**(3/4)*np.exp(-a*(x**2+y**2+z**2))
#     T=np.zeros((N, N))
#     for i in range(N):
#         for j in range(i, N):

# #コアハミルトニアン行列Hの計算
