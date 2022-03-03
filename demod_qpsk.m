function dmod = demod_qpsk(data)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
dmod = zeros(1,(2*length(data)));
demod_data = [];

for i = 1:length(data)
    
    
    data1 = norm(exp(j*pi/4) - data(i));
    data2 = norm(exp(j*3*pi/4) - data(i));
    data3 = norm(exp(j*5*pi/4) - data(i));
    data4 = norm(exp(j*7*pi/4) - data(i));
    
    tmp = [data1,data2,data3,data4];
    data_demod = min(tmp);
    d_demod = find(tmp == data_demod);
    
    [val, loc] = min(tmp);
    
    if (loc == 1)
        dmod(2*i) = 0;
        dmod(2*i -1) = 0;
    end
    if (loc == 2)
        
        dmod(2*i) = 0;
        dmod(2*i -1) = 1;
    end
    if (loc == 3)
        
        dmod(2*i) = 1;
        dmod(2*i -1) = 1;
    end
    if (loc == 4)
       
        dmod(2*i) = 1;
        dmod(2*i -1) = 0;
    end


end
return