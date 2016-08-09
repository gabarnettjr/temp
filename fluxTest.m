function fluxTest

stencilSize = 6;
phi = @(x,y)  (x.^2+y.^2) .^ (3/2);
poly = @(x,y) [ ones(length(x),1), x, y ];

n = 52;
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
    b = [ h^4*[2/3,1/6,2/3,2/3,1/6,2/3], [h,0,0] ];
  end
  w = b / A;
  W(:,i) = w(1:size(X,2));
end
ii = repmat( 1:size(X,1), size(X,2), 1 );
jj = idx.';
W = sparse( ii, jj, W, size(X,1), length(xe), size(X,2)*length(xe) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = poly( x, y )
z = [ ones(length(x),1), x, y ];
