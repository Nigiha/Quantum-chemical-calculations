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


        E_tot=RHF.RHF(temp_xyz, base_function_file, eps, max_iter, print_E=False)

        result.append({
            "variables": new_val,
            "energy": E_tot
        })

    return result



def optimize_H2O_type(raw_data, steps:list, val_ranges:list, base_function_file, eps, max_iter):

    a_locate_list=[]

    with open(raw_data, "r", encoding="utf-8") as f:
        lines=[line.split() for line in f]

    for i in range(2, len(lines)):
        for j in range(len(lines[i])):
            if lines[i][j]=="a":
                a_locate_list.append([i, j])
            
            elif lines[i][j]=="r":
                r_locate=[i, j]


    if len(a_locate_list)!=len(steps):
        raise ValueError("入力ファイル内の変数の数とstepsの数が一致しません")

    if len(a_locate_list)!=len(val_ranges):
        raise ValueError("入力ファイル内の変数の数とval_rangesの数が一致しません")
    
    loop=1
    grids=[]
    for (start, end), step in zip(val_ranges, steps):
        grid_array=np.arange(start, end+step, step)
        grids.append(grid_array)
        loop*=len(grid_array)


    grids_comb=list(itertools.product(*grids))

    result=[]

    loop_count=0

    min_E=float("inf")
    min_x=[]
    min_y=[]

    for new_val in grids_comb:#メインループ

        loop_count+=1
        print("ループ数: {:}/{:}".format(loop_count, loop))

        r=0

        for locate, new_data in zip(a_locate_list, new_val):
            i, j=locate
            lines[i][j]=str(new_data)

            r+=pow(new_data, 2)

        i, j=r_locate
        lines[i][j]=str(np.sqrt(r))

    
        temp_xyz="temp_scan.xyz"
        with open(temp_xyz, "w", encoding="utf-8") as f:
            for line in lines:
                f.write(" ".join(line)+"\n")


        E_tot=RHF.RHF(temp_xyz, base_function_file, eps, max_iter, print_E=True)
        
        if min_E>E_tot:
            min_E=E_tot
            min_x=[new_val[0]]
            min_y=[new_val[1]]
        elif min_E==E_tot:
            min_x.append(new_val[0])
            min_y.append(new_val[1])


        result.append({
            "variables": new_val,
            "energy": E_tot
        })

    return result, min_E, min_x, min_y


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

    ax.plot(a, E, color="red")

    plt.title("Structural_optimization of HeH+")
    plt.xlabel("The distance between He and H+ \n(Å)")
    plt.ylabel("Energy \n(hartree)")

    ax.grid()
    ax.legend()
    plt.show()

    

def double_d_figure_make(opt_result, min_E, min_x, min_y):
    plt.rcParams["backend"]="qtagg"
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1, projection="3d")

    x=[]
    y=[]
    E=[]

    for plot in opt_result:
        x.append(plot["variables"][0])
        y.append(plot["variables"][1])
        E.append(plot["energy"])

    surf=ax.plot_trisurf(x, y, E, cmap="viridis", edgecolor="none")
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=0.5, label="Energy (hartree)")

    for i in range(len(min_x)):
        x_val=min_x[i]
        y_val=min_y[i]

        r_val = np.sqrt(x_val**2 + y_val**2)
        theta_rad = np.arctan2(y_val, x_val)
        theta_dig=np.degrees(theta_rad)

        label_str = f"r={r_val:.3f} Å, θ={theta_dig:.1f}°\nE={min_E:.6f} hartree"

        ax.scatter(x_val, y_val, min_E, color="red", s=100, marker="*", label=label_str, zorder=5)

        

    z_bottom, z_top = ax.get_zlim()
    x_left, x_right = ax.get_xlim()
    ax.plot([0, 0], [0, 0], [z_bottom, z_top], color="black", linestyle="-.", linewidth=1.2, zorder=10)
    ax.text(0, 0, z_top, " r=0", color="black", fontsize=12, fontweight="bold", ha="left")
    ax.plot([0, x_right], [0, 0], [z_bottom, z_bottom], color="black", linestyle="-.", linewidth=1.2, zorder=10)
    ax.text(x_right, 0, z_bottom, " θ=0", color="black", fontsize=12, fontweight="bold", ha="left")

    ax.set_xlabel("x-coordinate (Å)")
    ax.set_ylabel("y-coordinate (Å)")
    ax.set_zlabel("Energy (hartree)")

    plt.title("Structural_optimization of H2O")
    # ax.set_zlabel("Energy \n(hartree)")

    ax.grid()
    ax.legend()
    plt.show()
    