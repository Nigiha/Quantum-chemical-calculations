import libs.optimize as opt
import pandas as pd
import time

start_time=time.time()

raw_data="samples/H2O_opt.xyz"
steps=[0.01, 0.01]
val_ranges=[[-0.40, 0.40], [0.80, 1.20]]
base_function_file="base_function_data/6-311G.json"
eps=1e-6
max_iter=100

opt_result, min_E, min_x, min_y=opt.optimize_H2O_type(raw_data, steps, val_ranges, base_function_file, eps, max_iter)

formatted_data = [
    {
        "x (Å)": res["variables"][0], 
        "y (Å)": res["variables"][1], 
        "Energy (hartree)": res["energy"]
    } 
    for res in opt_result
]
df=pd.DataFrame(formatted_data)
df.to_csv("Structural_optimization of H2O.csv", index=False)

end_time=time.time()
print("実行時間: {:.0f}s".format(end_time-start_time))

opt.double_d_figure_make(opt_result, min_E, min_x, min_y)
