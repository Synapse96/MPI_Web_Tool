function [compAO,compAOV1,compAOFinal]=getCompSigForAO(N,k2,w,k1,s1,s3,s5,s6)
% consturct composite signal for AO

compAO = zeros(1,N);
compAOV1 = zeros(1,N);
compAOFinal = zeros(1,N);

for i = (k2+w+1):(N-w-k2)
    
    id1 = (i-w-k2):(i-w-k1);
    id2 = (i+w+k1):(i+w+k2);    
    id3 = (i-w-k1):(i-w-1);
    id4 = (i+w+1):(i+w+k1);
    
    compAO(i) = s1(i)+2*mean(s1(id2))-2*mean(s3(id1))-mean(s3(id2));
    
end

compAO(compAO<0)=0;
compAO([1:k2+w N-w-k2+1:end])= 2*s1([1:k2+w N-w-k2+1:end]);

s3v2 = s3-s1;
s3v2(s3v2<0)=0;

for i = (k2+w+1):(N-w-k2)
    
    id1 = (i-w-k2):(i-w-k1);
    id2 = (i+w+k1):(i+w+k2);
    
    id3 = (i-w-k1):(i-w-1);
    id4 = (i+w+1):(i+w+k1);
        
    compAOV1(i)=compAO(i)-2*mean(compAO(id1))-2*mean(s3v2(id1))-mean(s3v2(id3));
end

compAOV1(compAOV1<0)=0;
compAOV1 = medfilt1(compAOV1,20);

compAOV2 = s5;
compAOV2(compAOV1==0)=0;


for i = (k2+w+1):(N-w-k2)
    
    id1 = (i-w-k2):(i-w-k1);
    id2 = (i+w+k1):(i+w+k2);   
    id3 = (i-w-k1):(i-w-1);
    id4 = (i+w+1):(i+w+k1);
    
    compAOFinal(i)=compAOV2(i)-mean(compAOV2(id3))-s6(i);
    
end

compAOFinal(compAOFinal<0)=0;