function result = calculate_a_interval(SignalFilePath)
	dataTemp = h5read(SignalFilePath, '/MovieGroup2/AcqPWCW/RawData/RawDataUnit');
    time = h5read(SignalFilePath, '/MovieGroup2/AcqPWCW/RawData/TimeStamp');
    
    dataNew = zeros(size(dataTemp));
    dataNew(1:64,:)=dataTemp(65:128,:);
    dataNew(65:128,:)=dataTemp(1:64,:);
    [addevent,aInterval] = RightMPI_a_wave(dataNew,time)
    result = aInterval(3,:)