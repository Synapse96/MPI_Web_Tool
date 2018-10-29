function ac =  findAC(s4,d3,roi)

x = filtfilt(d3,s4);
wid = x(roi);
[maxval,ac] = max(wid);
ac = roi(ac);




