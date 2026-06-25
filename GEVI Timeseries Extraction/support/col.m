
function c=col(num)
% DEFINE COLOR PROGRESSION
num = num+1;

c=[0 0 0];
num=mod(num,15);
switch num
    case 1
        c=[0.4 0.4 0.4];     
    case 2
        c=[0 160 0]/255;   
    case 3
        c=[0 0 255]/255;     
    case 4
        c=[104 34 139]/255;  
    case 5
        c=[255 0  0]/255;   
    case 6
        c=[225 118 0]/255;   
    case 7  
        c=[240 220 0]/255;  
    case 8
        c=[255 175 0]/255;   
    case 9
        c=[255 131 250]/255;
    case 10
        c=[205 0 205]/255;  
    case 11
        c=[50 153 204]/255;
    case 12
        c=[112 219 147]/255;     
    case 13
        c=[151 105 79]/255;   
    case 14
        c=[107 66 38]/255;   
    case 0
        c=[140 140 140]/255;  
end