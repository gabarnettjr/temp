function fluxTest

stencilSize = 6;
phi = @(x,y)  (x.^2+y.^2) .^ (3/2);
poly = @(x,y) [ ones(length(x),1), x, y ];

n = 52;                                                 % total number of cells going across the domain
h = 1/(n-2);                                            % width and height of one cell
[xx,yy] = meshgrid( -h/2 : h : 1+h/2 );                 % location of cell-averaged values in mesh form
x = xx(:);  y = yy(:);                                  % location of cell-averaged values in vector form
[xV,yV] = meshgrid( 0:h:1, h/2:h:1-h/2 );
xV = xV(:);  yV = yV(:);                                % location of vertical cell walls
[xH,yH] = meshgrid( h/2:h:1-h/2, 0:h:1 );
xH = xH(:);  yH = yH(:);                                % location of horizontal cell walls

idxV = knnsearch( [x,y], [xV,yV], 'k', stencilSize );   % index of nearest cell-average neighbors to each vertical wall
idxH = knnsearch( [x,y], [xH,yH], 'k', stencilSize );   % index of nearest cell-average neighbors to each horizontal wall

WV = quadMatrix( phi, poly, x, y, idxV, xV, yV, 1 );    % sparse matrix of quadrature weights along vertical cell walls
WH = quadMatrix( phi, poly, x, y, idxH, xH, yH, 0 );    % sparse matrix of quadrature weights along horizontal cell walls

indL = 1 : (n-2)^2;
indR = n-2+1 : (n-2)*(n-1);
indB = [];
indT = [];
for i = 1 : n-2
    indB = [ indB, (i-1)*(n-1)+(1:n-2) ];
    indT = [ indT, (i-1)*(n-1)+(2:n-1) ];
end

WVlr = WV(indL,:) - WV(indR,:);
WHbt = WH(indB,:) - WH(indT,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = quadMatrix( phi, poly, x, y, idx, xe, ye, vert )
X = x(idx);  Y = y(idx);
W = zeros( size(X,2), size(X,1) );
for i = 1 : size(X,1)
    xn = X(i,:) - xe(i);  yn = Y(i,:) - ye(i);
    xx = meshgrid(xn);  yy = meshgrid(yn);
    P = poly(xn,yn);
    A = [ phi(xx.'-xx,yy.'-yy), P;  P.', zeros(size(P,2),size(P,2)) ];
    if vert == 1
        b = [ h^4*[2/3,1/6,2/3,2/3,1/6,2/3], [h,0,0] ];
    else
        b = [ h^4*[2/3,2/3,1/6,1/6,2/3,2/3], [h,0,0] ];
    end
    w = b / A;
    W(:,i) = w( 1 : size(X,2) );
end
ii = repmat( 1:size(X,1), size(X,2), 1 );
jj = idx.';
W = sparse( ii, jj, W, size(X,1), length(xe), size(X,2)*length(xe) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = odefun( t, psi, u, v, x, y, WVlr, WHbt )
% z = WV(indL,:)*psiU - WV(indR,:)*psiU + WH(indB,:)*psiV - WH(indT,:)*psiV;
z = WVlr*(psi.*u(x,y,t)) + WHbt*(psi.*v(x,y,t));
