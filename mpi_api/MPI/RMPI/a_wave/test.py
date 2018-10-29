import matlab.engine
eng = matlab.engine.start_matlab()
SignalFilePath = 'C:/comp_files/Thesis_Project/IMG_20180413_16_64.raw'
print(eng.calculate_a_interval(SignalFilePath))