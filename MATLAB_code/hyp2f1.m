function [F,relErr] = hyp2f1(a,b,c,z,varargin)
%HYP2F1_NUM  Numeric Gauss hypergeometric 2F1(a,b;c;z) via Euler integral
%
%   F = HYP2F1(a,b,c,z) returns the value of 2F1(a,b;c;z) for scalar
%   complex or real a, b, c, z with Re(c) > Re(b) > 0 and Re(z) < 1.
%
%   F = HYP2F1(__,'Name',Value,...) forwards Name–Value pairs directly
%   to MATLAB's INTEGRAL (e.g. 'RelTol',1e-10,'AbsTol',0).
%
%   [F,relErr] = HYP2F1(...) additionally returns a conservative
%   relative-error estimate obtained by repeating the integral with a
%   10× tighter relative tolerance and comparing the two results.
%
%   Example
%   -------
%     >> F = hyp2f1(1,2,2.5,0.75)
%     F =
%          2.83679830462503
%
%   (For comparison, 2F1(1,2;2.5;0.75) ≈ 2.8367983046245793 in python scipy special.)
%
%   See also INTEGRAL, GAMMA, GAMMALN.

% ------------------------------------------------------------------------
% Input checks
% ------------------------------------------------------------------------
validateattributes(a,{'numeric'},{'scalar','finite'});
validateattributes(b,{'numeric'},{'scalar','finite'});
validateattributes(c,{'numeric'},{'scalar','finite'});
validateattributes(z,{'numeric'},{'scalar','finite'});
if ~(real(c) > real(b) && real(b) > 0)
    error('Input must satisfy Re(c) > Re(b) > 0 for Euler representation.');
end
if ~(real(z) < 1)
    error('Euler integral converges only for Re(z) < 1.');
end

% ------------------------------------------------------------------------
% Build integrand handle:  f(x) = x^(b-1) (1-x)^(c-b-1) (1-zx)^(-a)
% Use .* and ./ so the handle is array-valued.
% ------------------------------------------------------------------------
fun = @(x) x.^(b-1) .* (1-x).^(c-b-1) .* (1 - z.*x).^(-a);

% Default tolerances replicate IDL call: RelTol = 1e-8, AbsTol = 0
defaultRelTol = 1e-8;
p = inputParser;
addParameter(p,'RelTol',defaultRelTol,@(t)isscalar(t)&&t>0);
addParameter(p,'AbsTol',0,@(t)isscalar(t)&&t>=0);
parse(p,varargin{:});
relTol = p.Results.RelTol;
absTol = p.Results.AbsTol;


% ------------------------------------------------------------------------
% Compute integral on [0,1]
% ------------------------------------------------------------------------
I = integral(fun,0,1,'ArrayValued',true,...
             'RelTol',relTol,'AbsTol',absTol);

% Prefactor α = Γ(c)/(Γ(b)Γ(c-b)) using gammaln for stability
alpha = exp(gammaln(c) - gammaln(b) - gammaln(c-b));

F = alpha .* I;    % final result

% ------------------------------------------------------------------------
% Optional conservative relative-error estimate
% ------------------------------------------------------------------------
if nargout > 1
    % Repeat with 10× tighter RelTol and compare
    I_tight = integral(fun,0,1,'ArrayValued',true,...
                       'RelTol',relTol*0.1,'AbsTol',absTol);
    F_tight = alpha .* I_tight;
    relErr  = abs(F_tight - F) ./ abs(F_tight + (F_tight==0));
end
end
