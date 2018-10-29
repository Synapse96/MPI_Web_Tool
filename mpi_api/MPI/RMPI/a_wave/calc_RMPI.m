function rmpi = calc_RMPI(a,b)
	if isa(a,'double') & isa(b,'double')
		disp(a)
		disp(b)
		c = (a - b)/b
		rmpi = c
	end
