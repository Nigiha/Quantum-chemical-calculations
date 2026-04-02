import numpy as np
import scipy.special as sp

import libs.integrals as integrals

#==========step1:対象分子の設定==========


#======================================

total_e=10 #全電子数

target_molecule_file="samples/HeH+.xyz"
base_function_file="base_function_data/STO-3G.json"

#======================================


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

    print("対象の分子 : "+str(lines[1].strip()))

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
# def V_nn_matrix(molecule):

#     V=np.zeros((N, N))
    
#     for i in range(N):
#         for j in range(i+1, N):
#             r_i, r_j=np.array(molecule[i]["coord"]), np.array(molecule[j]["coord"])
#             Z_i, Z_j=symbol_to_Z[molecule[i]["symbol"]], symbol_to_Z[molecule[j]["symbol"]]
#             V[i][j]=Z_i*Z_j/(np.linalg.norm(r_i-r_j))
#             V[j][i]=Z_i*Z_j/(np.linalg.norm(r_i-r_j))
    
#     return V

def V_nn_schalar(molecule):

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
# 重なり行列Sの計算(GTOをs型に限れば、解析計算は簡単になる)

def S_prim(a:float, b:float, A:np.ndarray, B:np.ndarray): #原始GTO同士の重なり積分　←　多中心積分の公式とガウス積分から
        return (np.pi/(a+b))**(3/2) * np.exp(-a*b/(a+b)*np.linalg.norm(A-B)**2)

def S(molecule_basis):
    S_mat=np.zeros((K, K))

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
            
            S_mat[n][m]=s
            S_mat[m][n]=s
    
    return S_mat


#Sの対角化
val_S, vec_S=np.linalg.eigh(S(molecule_basis))

val_S_inv_sqrt = 1.0 / np.sqrt(val_S)
s_inv_sqrt=np.diag(val_S_inv_sqrt, 0)
U=vec_S


#変換行列Xの計算(正準直交化)                 
X=U@s_inv_sqrt


#==========step3:コアハミルトニアン行列の計算==========
#運動エネルギー項Tの計算
def T_prim(a:float, b:float, A:np.ndarray, B:np.ndarray):
    return a*b/(a+b)*(3-2*a*b/(a+b)*np.linalg.norm(A-B)**2)*S_prim(a, b, A, B)

def T(molecule_basis):
    T_mat=np.zeros((K, K))

    for m in range(K):
        for n in range(m, K):
            alpha_m, alpha_n=molecule_basis[m]["exponents"], molecule_basis[n]["exponents"]
            c_m, c_n=molecule_basis[m]["coefficients"], molecule_basis[n]["coefficients"]

            A_m=molecule_basis[m]["center"]
            A_n=molecule_basis[n]["center"]

            t=0

            for i in range(len(alpha_m)):
                for j in range(len(alpha_n)):
                        t+=c_m[i]*c_n[j]*Nm_c(alpha_m[i])*Nm_c(alpha_n[j])*T_prim(alpha_m[i], alpha_n[j], A_m, A_n)
            
            T_mat[m][n]=t
            T_mat[n][m]=t

    return T_mat

# print(T(molecule_basis))



#核引力項V_neの計算
def V_ne(molecule_basis):
    V_ne_mat=np.zeros((K, K))

    for m in range(K):
        for n in range(m, K):
            alpha_m, alpha_n=molecule_basis[m]["exponents"], molecule_basis[n]["exponents"]
            c_m, c_n=molecule_basis[m]["coefficients"], molecule_basis[n]["coefficients"]
            A_m, A_n=molecule_basis[m]["center"], molecule_basis[n]["center"]

            v=0

            for i in range(len(alpha_m)):
                for j in range(len(alpha_n)):
                    for atom in molecule:
                        atom_center=atom["coord"]
                        Z_c = symbol_to_Z[atom["symbol"]]
                        v+=(-Z_c)*c_m[i]*c_n[j]*Nm_c(alpha_m[i])*Nm_c(alpha_n[j])*integrals.nuclear_attraction(alpha_m[i], (0,0,0), A_m, alpha_n[j], (0,0,0), A_n, atom_center)
            
            V_ne_mat[m][n]=v
            V_ne_mat[n][m]=v

    return V_ne_mat



#コアハミルトニアン行列Hの計算
H=T(molecule_basis)+V_ne(molecule_basis)
# H=[[-2.62490, -1.50875],
#    [-1.50875, -1.77424]]

#==========step4:密度行列の初期値の設定==========

#直交化基底に対するコアハミルトニアン行列の計算
H_prime=X.T@H@X


#係数行列Cの計算(C_primeの初期値はH^を対角化する行列にしておく.最終的にはFock行列を対角化するものにしたい)
H_val, C_prime=np.linalg.eigh(H_prime)
C=X@C_prime

num_occupied=total_e//2 #占有軌道の数
C_occ=C[:, 0:num_occupied]


#密度行列Pの計算
P=2*C_occ@C_occ.T #P_(\mu\nu)=2\Sigma_i c_(\mu i)* c_(\nu i)






#==========step5:電子反発積分の計算==========
#(μν|λσ)の計算
#後回し
V_ee=[[1.0562, 0.4677, 0.6064],
      [0.4677, 0.2465, 0.3887],
      [0.6064, 0.3887, 0.7750]]

def calc_ERI(m, n, l, s):
    alpha_m, alpha_n, alpha_l, alpha_s=molecule_basis[m]["exponents"], molecule_basis[n]["exponents"], molecule_basis[l]["exponents"], molecule_basis[s]["exponents"]
    c_m, c_n, c_l, c_s=molecule_basis[m]["coefficients"], molecule_basis[n]["coefficients"], molecule_basis[l]["coefficients"], molecule_basis[s]["coefficients"]
    A_m, A_n, A_l, A_s=molecule_basis[m]["center"], molecule_basis[n]["center"], molecule_basis[l]["center"], molecule_basis[s]["center"]

    lmn=(0, 0, 0)

    eri=0

    for i_m in range(len(alpha_m)):
        for i_n in range(len(alpha_n)):
            for i_l in range(len(alpha_l)):
                for i_s in range(len(alpha_s)):
                    eri+=c_m[i_m]*c_n[i_n]*c_l[i_l]*c_s[i_s]*Nm_c(alpha_m[i_m])*Nm_c(alpha_n[i_n])*Nm_c(alpha_l[i_l])*Nm_c(alpha_s[i_s])*integrals.electron_repulsion(alpha_m[i_m], lmn, A_m, alpha_n[i_n], lmn, A_n, alpha_l[i_l], lmn, A_l, alpha_s[i_s], lmn, A_s)

    return eri






E_0_RHF_list=[float("inf")]
E_conv=1e-6
max_iter=100

for iteration in range(max_iter):
    #==========step6:Fock行列の計算==========

    #Coulomb積分Jの計算
    def coulomb_integral(m, n):
        I=0
        for l in range(K):
            for s in range(K):
                I+=calc_ERI(m, n, l, s)*P[l][s]
        return I

        
    #交換積分Kの計算
    def exchange_integral(m, n):
        I=0
        for l in range(K):
            for s in range(K):
                I+=calc_ERI(m, s, l, n)*P[l][s]
        return I


    #2電子項Gの計算
    G=np.zeros((K, K))
    for m in range(K):
        for n in range(K):
            G[m][n]=coulomb_integral(m, n)-(1/2)*exchange_integral(m, n)


    #Fock行列Fの計算
    F=H+G






    #===========step7:RHFエネルギーの計算==========

    #E_0^RHFの計算 (行列を用いた計算に書き換えられる)
    E_0_RHF=0
    for m in range(K):
        for n in range(K):
            E_0_RHF+=(1/2)*P[m][n]*(H[m][n]+F[m][n])

    E_0_RHF_list.append(E_0_RHF)

    #E_tot^RHFの計算
    E_tot_RHF=E_0_RHF+V_nn_schalar(molecule)


    #==========step8:収束判定==========
    #エネルギー変化による判定
    if abs(E_0_RHF_list[-1]-E_0_RHF_list[-2])<E_conv:
        break #step10へ



    #==========step9:Roothaan方程式の解法==========

    #直交化基底に対するFock行列F_primeの計算
    F_prime=X.T@F@X

    #C_primeの計算
    E, C_prime=np.linalg.eigh(F_prime)

    C=X@C_prime

    num_occupied=total_e//2 #占有軌道の数
    C_occ=C[:, 0:num_occupied]


    #密度行列Pの計算
    P=2*C_occ@C_occ.T #P_(\mu\nu)=2\Sigma_i c_(\mu i)* c_(\nu i)

else:
    print("not converge")



#==========step10:エネルギー値の出力==========
print("Total Energy : {:.6f} hartree".format(E_tot_RHF))
