function [ fVal ] = hFitFunc(a)
u = [5    10    15    20    25    30    35    40    45    50    55    60    65    70    75   100];
h=[
    0.1071    0.1087    0.2220    0.2419    0.2792    0.2313    0.3678    0.4205    0.3813    0.4521    0.4369    0.4481    0.4608    0.4762    0.4882    0.3952];
fVal = sqrt(sum((1-exp(-a(1)*u)./(1-a(2).*(exp(-a(1).*u)-1)) -h).^2));

end

