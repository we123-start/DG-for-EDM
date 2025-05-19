function U = minmod(a,b,c)
if sign(a)==sign(b)&&sign(b)==sign(c)
    d=min(abs(a),abs(b));
    U=sign(a)*min(d,abs(c));
else
    U=0;
end