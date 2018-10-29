import matlab.engine
eng = matlab.engine.start_matlab()
SignalFilePath = 'C:/comp_files/Thesis_Project/IMG_20180413_16_65.raw'
print(eng.calculate_mpi(SignalFilePath))