import numpy as np
import itertools
import matplotlib.pyplot as plt

import libs.RHF as RHF

#入力ファイル内の"a"のところに数値を代入して、各値でRHFを回す
def optimize(raw_data, steps:list, val_ranges:list, base_function_file, eps, max_iter):

    a_locate_list=[]

    with open(raw_data, "r", encoding="utf-8") as f:
        lines=[line.split() for line in f]

    for i in range(2, len(lines)):
        for j in range(len(lines[i])):
            if lines[i][j]=="a":
                a_locate_list.append([i, j])


    if len(a_locate_list)!=len(steps):
        raise ValueError("入力ファイル内の変数の数とstepsの数が一致しません")

    if len(a_locate_list)!=len(val_ranges):
        raise ValueError("入力ファイル内の変数の数とval_rangesの数が一致しません")
    
    
    grids=[]
    for (start, end), step in zip(val_ranges, steps):
        grids.append(np.arange(start, end+step, step))

    grids_comb=list(itertools.product(*grids))

    result=[]

    for new_val in grids_comb:
        for locate, new_data in zip(a_locate_list, new_val):
            i, j=locate
            lines[i][j]=str(new_data)

    
        temp_xyz="temp_scan.xyz"
        with open(temp_xyz, "w", encoding="utf-8") as f:
            for line in lines:
                f.write(" ".join(line)+"\n")


        E_tot=RHF.RHF(temp_xyz, base_function_file, eps, max_iter)

        result.append({
            "variables": new_val,
            "energy": E_tot
        })

    return result



#optimizeの結果の図を作成
def single_d_figure_make(opt_result):
    plt.rcParams["backend"]="qtagg"
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)

    a=[]
    E=[]

    for plot in opt_result:
        a.append(plot["variables"][0])
        E.append(plot["energy"])

    print(a, E)
    ax.plot(a, E, color="red")

    plt.title("Structural_optimization of HeH+")
    plt.xlabel("The distance between He and H")
    plt.ylabel("Energy")

    ax.grid()
    ax.legend()
    plt.show()

    


    




