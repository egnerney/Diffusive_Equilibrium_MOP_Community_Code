% quick_test_U.m
cases = [ ...
    1.5, 2.3, 4.0 ;   % moderate 
    0.75, 1.2, 0.8 ;  % small z
    3.0,  5.5, 7.0 ]; % larger a

for k = 1:size(cases,1)
    a = cases(k,1); b = cases(k,2); z = cases(k,3);
    val = hyperu(a,b,z,'RelTol',1e-10);
    fprintf('U(%g,%g,%g) = %.15g\n',a,b,z,val);
end