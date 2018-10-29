function result = calculate_mpi(SignalFilePath)
	dataTemp = h5read(SignalFilePath, '/MovieGroup2/AcqPWCW/RawData/RawDataUnit');
    time = h5read(SignalFilePath, '/MovieGroup2/AcqPWCW/RawData/TimeStamp');
    
    dataNew = zeros(size(dataTemp));
    dataNew(1:64,:)=dataTemp(65:128,:);
    dataNew(65:128,:)=dataTemp(1:64,:);
    [resultTable,resultValTime,MCnew,valueMC,AOnew,valueAO,ACnew,valueAC,MOnew,valueMO] = mainfile1028F(dataNew,time)
    result = resultTable(:,8)