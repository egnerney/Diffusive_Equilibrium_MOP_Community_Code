% quick_test.m   — simple sanity-check for hyp2f1.m

vec = [ ...
    1   , 2   , 2.5 , 0.75 ;   % reference value ≈ 2.8367983046245793 (scipy special)
    0.5 , 1   , 1.8 , 0.50 ;   % any valid combo (Re(c)>Re(b)>0, Re(z)<1)
    2.5 , 3.5 , 5   , 0.20 ];

for k = 1:size(vec,1)
    row = vec(k,:);
    a = row(1); b = row(2); c = row(3); z = row(4);

    val = hyp2f1(a,b,c,z,'RelTol',1e-10);   % call our function
    fprintf('2F1(%g,%g;%g;%g) = %.15g\n',a,b,c,z,val);
end