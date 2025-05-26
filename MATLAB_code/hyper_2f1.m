function fz = hyper_2f1 ( a, b, c, z )

%*****************************************************************************80
%
%% hyper_2f1() evaluates the hypergeometric function 2F1(a,b,c,z).
%
%  Discussion:
%
%    This is simply a convenient interface to the built-in hypergeom() function.
%
%  Licensing:
%
%    This code is distributed under the MIT license.
%
%  Modified:
%
%    22 December 2023
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    real a, b, c: the parameters.
%
%    real or complex z: the argument.
%
%  Output:
%
%    real or complex fz: the function value.
%
  fz = hypergeom ( [ a, b ], c, z );

  return
end

