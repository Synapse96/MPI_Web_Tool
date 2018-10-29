function result = pad_arr(arr,len)
	[m,n] = size(arr) 
	arr(n,n + len) = {0}
	result = arr