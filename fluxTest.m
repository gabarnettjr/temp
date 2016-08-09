function fluxTest

stencilSize = 6;
phi = @(x,y)  (x.^2+y.^2) .^ (3/2);

n = 5;
h = 1/(n-2);
[xx,yy] = meshgrid( -h/2 : h : 1+h/2 );
x = xx(:);  y = yy(:);
[xV,yV] = meshgrid( 0:h:1, h/2:h:1-h/2 );
xV = xV(:);  yV = yV(:);
[xH,yH] = meshgrid( h/2:h:1-h/2, 0:h:1 );
xH = xH(:);  yH = yH(:);

idxV = knnsearch( [x,y], [xV,yV], 'k', stencilSize );
idxH = knnsearch( [x,y], [xH,yH], 'k', stencilSize );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = quadMatrixV( phi, x, y, idxV, xV, yV )
W = zeros( size(idxV,2), length(xV) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = quadMatrixH( phi, x, y, idxH, xH, yH )
W = zeros( size(idxH,2), length(xH) );
