function out = linear( x )
if x<0.5
    out = 2.*x;
else 
    out = 2.*(1-x);
end