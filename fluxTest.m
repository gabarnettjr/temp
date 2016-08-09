% Conservative transport scheme using PHS RBFs and polynomials for integration across cell walls.
% The test problem is the deforming bubble from Blossey and Durran (2008)
% Greg Barnett
% August 9th, 2016

function fluxTest

stencilSize = 6;                                        % number of nodes to use per flux integral calculation (don't change this yet)
phi = @(x,y)  (x.^2+y.^2) .^ (3/2);                     % RBF which will be shifted to create RBF part of basis
poly = @(x,y) [ ones(1,length(x)); x; y ];              % polynomial basis functions
n = 52;                                                 % total number of cells going across the domain, including 2 ghost cells
t = 0 : 1/200 : 1;                                      % vector of time values

h = 1/(n-2);                                            % width and height of one cell (space step)
k = t(2) - t(1);                                        % time elapsed during one time step
[xx,yy] = meshgrid( -h/2 : h : 1+h/2 );
rtilde = 5 * sqrt( (xx-.3).^2 + (yy-.5).^2 );
psi = ( (1+cos(pi*rtilde)) ./ 2 ) .^ 2;                 % initial condition for scalar field psi
psi(rtilde>1) = 0;

x = xx(:);  y = yy(:);                                  % location of cell-averaged values
[xV,yV] = meshgrid( 0:h:1, h/2:h:1-h/2 );
xV = xV(:);  yV = yV(:);                                % location of vertical cell walls (midpoint)
[xH,yH] = meshgrid( h/2:h:1-h/2, 0:h:1 );
xH = xH(:);  yH = yH(:);                                % location of horizontal cell walls (midpoint)

idxV = knnsearch( [x,y], [xV,yV], 'k', stencilSize );   % index of nearest cell-average neighbors to each vertical wall
idxH = knnsearch( [x,y], [xH,yH], 'k', stencilSize );   % index of nearest cell-average neighbors to each horizontal wall

WV = quadMatrix( phi, poly, h, x, y, idxV, xV, yV, 1 ); % sparse matrix of quadrature weights along vertical cell walls
WH = quadMatrix( phi, poly, h, x, y, idxH, xH, yH, 0 ); % sparse matrix of quadrature weights along horizontal cell walls

indL = 1 : (n-2)^2;                                     % index for left wall of each interior cell
indR = n-2+1 : (n-2)*(n-1);                             % index for right wall of each interior cell
indB = [];
indT = [];
for i = 1 : n-2
    indB = [ indB, (i-1)*(n-1)+(1:n-2) ];               % index for bottom wall of each interior cell
    indT = [ indT, (i-1)*(n-1)+(2:n-1) ];               % index for top wall of each interior cell
end

WVlr = WV(indL,:) - WV(indR,:);                         % combine to get total flux through left and right walls
WHbt = WH(indB,:) - WH(indT,:);                         % combine to get total flux through bottom and top walls

U = @(t) u(x,y,t);  V = @(t) v(x,y,t);
f = @(t,psi)  odefun( t, psi, U(t), V(t), WVlr, WHbt );

for i = 1 : length(t)-1
    psi = rk( t(i), psi, k, f, 4 );
    if mod( t(i)*100, 5) == 0
        psi = reshape( psi, n, n );
        figure(1),clf
            contour( xx, yy, psi, -.05:.1:.95 )
            axis( 'equal', [0,1,0,1] )
            caxis( [-.1,1] )
            colorbar( parula(10) )
        drawnow
        psi = psi(:);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = quadMatrix( phi, poly, h, x, y, idx, xe, ye, vert )
stencilSize = size(idx,2);
ne = length(xe);
X = x(idx);  Y = y(idx);
W = zeros( stencilSize, ne );
for i = 1 : ne
    xn = X(i,:) - xe(i);  yn = Y(i,:) - ye(i);
    xx = meshgrid(xn);  yy = meshgrid(yn);
    P = poly(xn,yn);
    A = [ phi(xx.'-xx,yy.'-yy), P.';  P, zeros(size(P,1),size(P,1)) ];
    if vert == 1
        b = [ h^4*[2/3,1/6,2/3,2/3,1/6,2/3], [h,0,0] ];
    else
        b = [ h^4*[2/3,2/3,1/6,1/6,2/3,2/3], [h,0,0] ];
    end
    w = b / A;
    W(:,i) = w( 1 : stencilSize );
end
ii = repmat( 1:ne, stencilSize, 1 );
W = sparse( ii, idx.', W, ne, length(x), stencilSize*ne );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% U and V should be functions of t only.

function z = odefun( t, psi, U, V, WVlr, WHbt )
psi = reshape( psi, sqrt(length(psi)), sqrt(length(psi)) );
psi([1,end],:) = psi([2,end-1],:);
psi(:,[1,end]) = psi(:,[2,end-1]);
psi = psi(:);
% z = WV(indL,:)*psiU - WV(indR,:)*psiU + WH(indB,:)*psiV - WH(indT,:)*psiV;
z = WVlr*(psi.*U(t)) + WHbt*(psi.*V(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = u(x,y,t)
T = 1;
r = sqrt( (x-.5).^2 + (y-.5).^2 );
th = atan2( y-.5, x-.5 );
uth = 4*pi*r./T .* ( 1 - cos(2*pi*t./T) .* (1-(4*r).^6)./(1+(4*r).^6) );
z = uth .* sin(th);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = v(x,y,t)
T = 1;
r = sqrt( (x-.5).^2 + (y-.5).^2 );
th = atan2( y-.5, x-.5 );
uth = 4*pi*r./T .* ( 1 - cos(2*pi*t./T) .* (1-(4*r).^6)./(1+(4*r).^6) );
z = -uth .* cos(th);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



























