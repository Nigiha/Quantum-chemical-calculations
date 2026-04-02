import numpy as np

import libs.integrals as integrals

#核間反発エネルギー
def V_nn_schalar(molecule_data):

    N=molecule_data["N"]
    molecule=molecule_data["atom_shells"]

    V=0
    
    for i in range(N):
        for j in range(i+1, N):
            r_i, r_j=np.array(molecule[i]["coord"]), np.array(molecule[j]["coord"])
            Z_i, Z_j=molecule[i]["Z"], molecule[j]["Z"]
            V+=Z_i*Z_j/(np.linalg.norm(r_i-r_j))
    
    return V


# def V_nn_mat(molecule_data):

#     N=molecule_data["N"]
#     molecule=molecule_data["atom_shells"]

#     V=np.zeros((N, N))
   
#     for i in range(N):
#         for j in range(i+1, N):
#             r_i, r_j=np.array(molecule[i]["coord"]), np.array(molecule[j]["coord"])
#             Z_i, Z_j=molecule[i]["Z"], molecule[j]["Z"]
#             V[i][j]=Z_i*Z_j/(np.linalg.norm(r_i-r_j))
#             V[j][i]=Z_i*Z_j/(np.linalg.norm(r_i-r_j))
    
#     return V



def extract_shell_params(shell):
    return (
        shell["exponents"],
        shell["coefficients"],
        shell["normalization_factors"],
        shell["angular_momentums"],
        shell["center"]
    )



#重なり積分
def S_mat(molecule_basis):

    K=molecule_basis["K"]
    e_shells=molecule_basis["electron_shells"]

    S_mat=np.zeros((K, K))

    for m in range(K):
        for n in range(m, K):
            alpha_m, c_m, Nm_m, l_m, A_m = extract_shell_params(e_shells[m])
            alpha_n, c_n, Nm_n, l_n, A_n = extract_shell_params(e_shells[n])

            s=0

            for i in range(len(alpha_m)):
                for j in range(len(alpha_n)):
                    s+=c_m[i]*c_n[j]*Nm_m[i]*Nm_n[j]*integrals.overlap(alpha_m[i], l_m, A_m, alpha_n[j], l_n, A_n)
            
            S_mat[m][n]=s
            S_mat[n][m]=s
    
    return S_mat



#運動エネルギー項
def T_mat(molecule_basis):

    K=molecule_basis["K"]
    e_shells=molecule_basis["electron_shells"]

    T_mat=np.zeros((K, K))

    for m in range(K):
        for n in range(m, K):
            alpha_m, c_m, Nm_m, l_m, A_m = extract_shell_params(e_shells[m])
            alpha_n, c_n, Nm_n, l_n, A_n = extract_shell_params(e_shells[n])

            t=0

            for i in range(len(alpha_m)):
                for j in range(len(alpha_n)):
                        t+=c_m[i]*c_n[j]*Nm_m[i]*Nm_n[j]*integrals.kinetic(alpha_m[i], l_m, A_m, alpha_n[j], l_n, A_n)
            
            T_mat[m][n]=t
            T_mat[n][m]=t

    return T_mat



#核引力項
def V_ne_mat(molecule_basis, molecule_data):

    K=molecule_basis["K"]
    e_shells=molecule_basis["electron_shells"]
    molecule=molecule_data["atom_shells"]

    V_ne_mat=np.zeros((K, K))

    for m in range(K):
        for n in range(m, K):
            alpha_m, c_m, Nm_m, l_m, A_m = extract_shell_params(e_shells[m])
            alpha_n, c_n, Nm_n, l_n, A_n = extract_shell_params(e_shells[n])

            v=0

            for i in range(len(alpha_m)):
                for j in range(len(alpha_n)):
                    for atom in molecule:
                        atom_center=atom["coord"]
                        Z_c = atom["Z"]
                        v+=(-Z_c)*c_m[i]*c_n[j]*Nm_m[i]*Nm_n[j]*integrals.nuclear_attraction(alpha_m[i], l_m, A_m, alpha_n[j], l_n, A_n, atom_center)
            
            V_ne_mat[m][n]=v
            V_ne_mat[n][m]=v

    return V_ne_mat



#電子反発項
def V_ee_tensor(molecule_basis):

    K=molecule_basis["K"]
    e_shells=molecule_basis["electron_shells"]

    V_ee_tens=np.zeros((K, K, K, K))

    for m in range(K):
        for n in range(m+1):

            mn=m*(m+1)//2+n #複合インデックス

            alpha_m, c_m, Nm_m, l_m, A_m = extract_shell_params(e_shells[m])
            alpha_n, c_n, Nm_n, l_n, A_n = extract_shell_params(e_shells[n])

            for l in range(K):
                for s in range(l+1):

                    ls=l*(l+1)//2+s

                    if mn < ls:
                        continue
                    
                    alpha_l, c_l, Nm_l, l_l, A_l = extract_shell_params(e_shells[l])
                    alpha_s, c_s, Nm_s, l_s, A_s = extract_shell_params(e_shells[s])

                    v=0

                    for i_m in range(len(alpha_m)):
                        for i_n in range(len(alpha_n)):
                            for i_l in range(len(alpha_l)):
                                for i_s in range(len(alpha_s)):
                                    v+=c_m[i_m]*c_n[i_n]*c_l[i_l]*c_s[i_s]* \
                                        Nm_m[i_m]*Nm_n[i_n]*Nm_l[i_l]*Nm_s[i_s]* \
                                        integrals.electron_repulsion(alpha_m[i_m], l_m, A_m, alpha_n[i_n], l_n, A_n, alpha_l[i_l], l_l, A_l, alpha_s[i_s], l_s, A_s)
                                    
                    V_ee_tens[m, n, l, s]=v
                    V_ee_tens[n, m, l, s]=v
                    V_ee_tens[m, n, s, l]=v
                    V_ee_tens[n, m, s, l]=v
                    V_ee_tens[l, s, m, n]=v
                    V_ee_tens[s, l, m, n]=v
                    V_ee_tens[l, s, n, m]=v
                    V_ee_tens[s, l, n, m]=v
            
    return V_ee_tens

#Coulomb積分Jの計算
def coulomb_integral(m, n, V_ee, P, K):
    I=0
    
    for l in range(K):
        for s in range(K):
            I+=V_ee[m, n, l, s]*P[l][s]
    
    return I

#交換積分Kの計算
def exchange_integral(m, n, V_ee, P, K):
    I=0

    for l in range(K):
        for s in range(K):
            I+=V_ee[m, s, l, n]*P[l][s]
    return I

#2電子項Gの計算
def G(V_ee, P, molecule_basis):
    K=molecule_basis["K"]
    G_mat=np.zeros((K, K))
    for m in range(K):
        for n in range(K):
            G_mat[m][n]=coulomb_integral(m, n, V_ee, P, K)-(1/2)*exchange_integral(m, n, V_ee, P, K)
    
    return G_mat