# import numpy as np
# from scipy import linalg
# from pyscf import gto

# mol=gto.M(
#     atom='H 0 0 0; F 0 0 0.91',
#     basis='sto-3g',
#     charge=0,
#     spin=0
# )

# nocc=mol.nelectron // 2
# Enuc=mol.energy_nuc()

# S=mol.intor('int1e_ovlp')
# T=mol.intor('int1e_kin')
# V=mol.intor('int1e_nuc')
# eri=mol.intor('int2e')
# Hcore=T+V

# evals, evecs=linalg.eigh(Hcore, S)
# C_occ=evecs[:, :nocc]

# P=2*np.einsum('ui, vi->uv', C_occ, C_occ)

# max_iter=50
# tol=1e-8
# old_energy=0.0

# print("Iter | Total Energy (hartree)  | dE")
# print("-"*45)

# for i in range(max_iter):
#     J=np.einsum('pqrs, rs->pq', eri, P)
#     K=np.einsum('prqs, rs->pq', eri, P)

#     F=Hcore+J-0.5*K

#     E_elec=0.5*np.sum(P*(Hcore+F))
#     E_tot=E_elec+Enuc

#     dE=abs(E_tot-old_energy)
#     print(f"{i:4d} | {E_tot:21.10f} | {dE:.3e}")

#     if dE<tol:
#         print("-"*45)
#         print("SCF Converged!")
#         break

#     old_energy=E_tot

#     evals, evecs=linalg.eigh(F, S)
#     C_occ=evecs[:, :nocc]
#     P=2*np.einsum('ui, vi->uv', C_occ, C_occ)

# else:
#     print("SCF did not converge.")


import numpy as np
from scipy import linalg
from pyscf import gto, dft

# 1. 分子のセットアップ
mol = gto.M(atom='H 0 0 0; F 0 0 0.91', basis='sto-3g', verbose=0)
nocc = mol.nelectron // 2
Enuc = mol.energy_nuc()

# --- DFT用の追加セットアップ ---
# グリッドの生成（動径方向・角度方向に点をばらまく）
grids = dft.gen_grid.Grids(mol)
grids.build()
# 各グリッド点における基底関数(AO)の値を事前計算 (点数 × 基底関数数 の行列)
ao_value = dft.numint.eval_ao(mol, grids.coords)
# -----------------------------

# 積分の取得 (1電子積分 + 2電子積分)
S = mol.intor('int1e_ovlp')
T = mol.intor('int1e_kin')
V = mol.intor('int1e_nuc')
eri = mol.intor('int2e')
Hcore = T + V

# 初期推測
evals, evecs = linalg.eigh(Hcore, S)
P = 2 * np.einsum('ui,vi->uv', evecs[:, :nocc], evecs[:, :nocc])

max_iter = 50
tol = 1e-8
old_energy = 0.0

print("Iter | Total Energy (Hartree) | dE")
print("-" * 45)

for i in range(max_iter):
    # クーロン行列 J の構築 (これはHFと同じ)
    J = np.einsum('pqrs,rs->pq', eri, P)
    
    # --- XC行列 (Vxc) の構築 ---
    # a) 各グリッド点での電子密度 rho を計算: rho = Σ_uv P_uv * phi_u * phi_v
    rho = np.einsum('pi,ij,pj->p', ao_value, P, ao_value)
    
    # b) 汎関数の評価 (ここではLDA汎関数の定番 'LDA,VWN' を使用)
    # exc: 1電子あたりのXCエネルギー, vxc: XCポテンシャル (vrho, vsigma, vlapl, vtau)
    exc, vxc_tuple, _, _ = dft.libxc.eval_xc('LDA,VWN', rho)
    vrho = vxc_tuple[0] # LDAでは密度微分(vrho)のみを使う
    
    # c) 数値積分による Vxc 行列の構築
    # Vxc_uv = Σ_p (重み_p * vrho_p) * phi_u(p) * phi_v(p)
    weight_vrho = grids.weights * vrho
    Vxc = np.einsum('pi,p,pj->ij', ao_value, weight_vrho, ao_value)
    
    # d) XCエネルギーの計算: Exc = ∫ rho(r) * exc(r) dr
    E_xc = np.sum(grids.weights * rho * exc)
    # ---------------------------
    
    # Kohn-Sham 行列 (Fock行列のDFT版) の構築
    F = Hcore + J + Vxc
    
    # 電子エネルギーと全エネルギーの計算
    # DFTのエネルギーはHFとは少し数式が異なります
    E_elec = np.sum(P * Hcore) + 0.5 * np.sum(P * J) + E_xc
    E_tot = E_elec + Enuc
    
    # 収束判定
    dE = abs(E_tot - old_energy)
    print(f"{i:4d} | {E_tot:21.10f} | {dE:.3e}")
    if dE < tol:
        print("-" * 45)
        print("DFT (LDA) Converged!")
        break
    old_energy = E_tot
    
    # 対角化と密度行列の更新
    evals, evecs = linalg.eigh(F, S)
    P = 2 * np.einsum('ui,vi->uv', evecs[:, :nocc], evecs[:, :nocc])

else:
    print("SCF did not converge.")