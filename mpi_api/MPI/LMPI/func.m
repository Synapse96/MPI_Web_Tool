function a = func(data)
	[Gx1, Gy1] = imgradientxy(imadjust(data));
	a = sum(Gx1(1:64,:))