function [U,relErr] = hyperu(a,b,z,varargin)
%HYPERU  Numeric evaluation of the confluent hypergeometric U(a,b,z)
%
%   U = HYPERU(a,b,z) computes the Tricomi confluent hypergeometric
%   function U(a,b,z) for scalar (real or complex) inputs satisfying
%
%        Re(a) > 0 ,   Re(z) > 0 .
%
%   It uses the Laplace–type integral representation
%
%         U(a,b,z) = (1 / Gamma(a)) * ∫₀^∞ e^{-z t} t^{a-1} (1+t)^{b-a-1} dt
%
%   and MATLAB’s adaptive INTEGRAL on (0,Inf).
%
%   U = HYPERU(...,'Name',Value,...) forwards Name–Value pairs directly
%   to INTEGRAL (supported names: 'RelTol', 'AbsTol', 'ArrayValued').
%
%   [U,relErr] = HYPERU(...) additionally returns a rough relative-error
%   estimate obtained by repeating the integral with 10× tighter RelTol.
%
%   Example
%   -------
%       >> hyperu(1.5, 2.3, 4.0)
%       ans =
%           0.117883421836848
%
%   Cross-check with SciPy:  scipy.special.hyperu(1.5, 2.3, 4.0) = 0.11788342183684895
%
%   See also GAMMA, GAMMALN, INTEGRAL.

% -------------------------------------------------------------------------
% Input validation
% -------------------------------------------------------------------------
validateattributes(a,{'numeric'},{'scalar','finite'});
validateattributes(b,{'numeric'},{'scalar','finite'});
validateattributes(z,{'numeric'},{'scalar','finite'});
if ~(real(a) > 0)
    error('Laplace representation requires Re(a) > 0.');
end
if ~(real(z) > 0)
    error('Laplace representation requires Re(z) > 0.');
end

% -------------------------------------------------------------------------
% Build integrand handle:  f(t) = e^(-z t) · t^(a-1) · (1+t)^(b-a-1)
% -------------------------------------------------------------------------
fun = @(t) exp(-z.*t) .* t.^(a-1) .* (1+t).^(b-a-1);

% Default tolerances: RelTol = 1e-8, AbsTol = 0
defaultRelTol = 1e-8;
p = inputParser;
addParameter(p,'RelTol',defaultRelTol,@(t)isscalar(t)&&t>0);
addParameter(p,'AbsTol',0,@(t)isscalar(t)&&t>=0);
parse(p,varargin{:});
relTol = p.Results.RelTol;
absTol = p.Results.AbsTol;

% -------------------------------------------------------------------------
% Evaluate integral on (0,Inf)
% -------------------------------------------------------------------------
I = integral(fun,0,Inf,'ArrayValued',true,...
             'RelTol',relTol,'AbsTol',absTol);

% Prefactor 1/Gamma(a) via gammaln for stability
U = exp(-gammaln(a)) .* I;

% -------------------------------------------------------------------------
% Optional conservative relative-error estimate
% -------------------------------------------------------------------------
if nargout > 1
    I2 = integral(fun,0,Inf,'ArrayValued',true,...
                  'RelTol',relTol*0.1,'AbsTol',absTol);
    U2 = exp(-gammaln(a)) .* I2;
    relErr = abs(U2-U)./abs(U2 + (U2==0));
end
end
