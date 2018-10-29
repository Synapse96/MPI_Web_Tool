import matlab.engine
import math 
eng = matlab.engine.start_matlab()
SignalFilePath = 'C:/comp_files/Thesis_Project/IMG_20180323_1_3.raw'
cell_array = eng.calculate_mpi(SignalFilePath)
res_dict = {'MPI': []}
res_dict['MPI'] = list(filter(lambda x: not math.isnan(x), cell_array[1:]))
print(res_dict)