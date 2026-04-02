import numpy as np

import libs.build_data as build_data
import libs.integrals as integrals
import libs.mtr_make as mtr



def RHF(input_file, base_function_file, eps, max_iter):

    #==========step1:対象分子の設定==========
    molecule_data=build_data.load_xyz(input_file, to_bohr=True)
    total_e=molecule_data["N_e"]
    num_occupied=total_e//2 #占有軌道の数

    molecule_basis=build_data.build_molecule_basis(molecule_data, base_function_file)
    K=molecule_basis["K"]

    V_nn=mtr.V_nn_schalar(molecule_data)

    #==========step2:正準直交化の実行==========
    S=mtr.S_mat(molecule_basis)
    #Sの対角化
    val_S, vec_S=np.linalg.eigh(S)

    val_S_inv_sqrt = 1.0 / np.sqrt(val_S)
    s_inv_sqrt=np.diag(val_S_inv_sqrt, 0)
    U=vec_S

    #変換行列Xの計算(正準直交化)                 
    X=U@s_inv_sqrt


    #==========step3:コアハミルトニアン行列の計算==========
    T=mtr.T_mat(molecule_basis)
    V_ne=mtr.V_ne_mat(molecule_basis, molecule_data)
    #コアハミルトニアン行列Hの計算
    H=T+V_ne


    #==========step4:密度行列の初期値の設定==========
    #直交化基底に対するコアハミルトニアン行列の計算
    H_prime=X.T@H@X

    #係数行列Cの計算(C_primeの初期値はH^を対角化する行列にしておく.最終的にはFock行列を対角化するものにしたい)
    H_val, C_prime=np.linalg.eigh(H_prime)
    C=X@C_prime

    C_occ=C[:, 0:num_occupied]

    #密度行列Pの計算
    P=2*C_occ@C_occ.T #P_(\mu\nu)=2\Sigma_i c_(\mu i)* c_(\nu i)


    #==========step5:電子反発積分の計算==========
    E_0_RHF_list=[float("inf")]

    V_ee=mtr.V_ee_tensor(molecule_basis)

    for iteration in range(max_iter):
        #==========step6:Fock行列の計算==========
        #2電子項Gの計算
        G=mtr.G(V_ee, P, molecule_basis)

        #Fock行列の計算
        F=H+G


        #===========step7:RHFエネルギーの計算==========
        #E_0^RHFの計算 
        E_0_RHF = np.sum(0.5 * P * (H + F))
        E_0_RHF_list.append(E_0_RHF)

        #E_tot^RHFの計算
        E_tot_RHF=E_0_RHF+V_nn


        #==========step8:収束判定==========
        #エネルギー変化による判定
        if abs(E_0_RHF_list[-1]-E_0_RHF_list[-2])<eps:
            break #step10へ

        #==========step9:Roothaan方程式の解法==========
        #直交化基底に対するFock行列F_primeの計算
        F_prime=X.T@F@X

        #C_primeの計算
        E, C_prime=np.linalg.eigh(F_prime)

        C=X@C_prime

        C_occ=C[:, 0:num_occupied]


        #密度行列Pの計算
        P=2*C_occ@C_occ.T #P_(\mu\nu)=2\Sigma_i c_(\mu i)* c_(\nu i)

    else:
        print("not converge")

    
    #==========step10:エネルギー値の出力==========
    print("Total Energy : {:.6f} hartree".format(E_tot_RHF))

    return E_tot_RHF
