# import libs.RHF as RHF
# input_file="samples/H2O.xyz"
# base_function_file="base_function_data/STO-3G.json"
# eps=1e-8
# max_iter=100
# RHF.RHF(input_file, base_function_file, eps, max_iter)





import libs.optimize as opt
raw_data="samples/HeH+_opt.xyz"
steps=[0.01]
val_ranges=[[0.5, 1.5]]
base_function_file="base_function_data/STO-3G.json"
eps=1e-6
max_iter=100

opt_result=opt.optimize(raw_data, steps, val_ranges, base_function_file, eps, max_iter)

opt.single_d_figure_make(opt_result)