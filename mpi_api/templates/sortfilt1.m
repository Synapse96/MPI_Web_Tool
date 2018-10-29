%y=sortfilt(x,n,p);
%applys a sorting filter to vector x of length n extracting the pth
%percentile element for y[n]
%
%Examples:
%minimum filter: y=sortfilt1(x,50,0);
%maximum filter: y=sortfilt1(x,50,100);
%median filter:  y=sortfilt1(x,50,50);
%'95th percentile' filter: y=sortfilt1(x,50,95);
%'10th percentile' filter: y=sortfilt1(x,50,10);
%Note if p<0 then p set to 0; if p>100 then p set to 100

%Philip de Chazal
%10th March 2006
%% this function take as imput a X vector of rowdata, n= lengh of the windows for computation, p= percentile.
%give as output a vector of the same lenght of X and whit the value of the percentile for evry samples computed with a windows of n length
% if p=0 we can see the vector of min val
% if p=100 we can see the vector of max val
% if p=50 we can see the vector of medium val computed in the window-n

%%

function y = sortfilt1(x,n,p)

%applys a sorting filter to vector x of window length n extracting the pth
%percentile element for y
%
%Examples:
%minimum filter: y=sortfilt1(x,50,0);
%maximum filter: y=sortfilt1(x,50,100);
%median filter:  y=sortfilt1(x,50,50);
%'95th percentile' filter: y=sortfilt1(x,50,95);
%'10th percentile' filter: y=sortfilt1(x,50,10);
%Note if p<0 then p set to 0; if p>100 then p set to 100

% Must be an integer window width
if n<0 || mod(n,1)~=0
    error('n must be an integer greater than 0');
end

% length of signal
N = length(x);

% ensure percentile value makes sense
if p>100
    p=100;
elseif p<0
    p=0;
end


% split window to center on current sample
% Different for odd of even window lengths
nEven = (mod(n,2)==0); % here we verify if n is a pair or odd number
if nEven
    N1 = (n/2)-1;
    N2 = (n/2);
else
    N1 = (n-1)/2;
    N2 = (n-1)/2;
end

% Set up output array
y=zeros(size(x));

for i=1:N
    
    range = max([1 i-N1]):min([N,i+N2]); % in this line we make the range whit a particular attention at the border condition!!!
    % ES. y = sortfilt1(x,n,p)
    % X=[1 2 3 4 5 6 7 8 9 10]
    % n=4 and p=100
    %---we can compute: N1=1 and N2=2 and N=10;
    % i=1...range=[max(1 1-1):min(10 1+2)]---range=[1:3];----y=3
    % i=2...range=[1:4];----y=4
    % i=3...range=[2:5];----y=5
    %[...]
    % i=9...range=[8:10];---y=10;
    % i=10..range=[9:10];---y=10;
    % y=[3 4 5 6 7 8 9 10 10 10];---we can see that the lenght of y is 10
    % equal to X lenght.
    
    % Round of to find which index in window is closest to percentile
    P = 1+ round((p/100)*(length(range)-1)); %Index to give particular percentile
%     disp([max([1 i-N1]) min([N,i+N2]) length(range) P])

    % select only those values in this range from sorted array
    Z = sort(x(range));
    
    y(i) = Z(P);
    
end