import RHF
input_file="samples/H2O.xyz"
base_function_file="STO-3G.json"
eps=1e-8
max_iter=100
RHF.RHF(input_file, base_function_file, eps, max_iter)