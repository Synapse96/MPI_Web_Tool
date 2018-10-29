import matlab.engine
eng = matlab.engine.start_matlab()
SignalFilePath = 'C:/comp_files/Thesis_Project/IMG_20180323_1_2.raw'
print(eng.calculate_mpi(SignalFilePath))