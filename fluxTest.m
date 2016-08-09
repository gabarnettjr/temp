clc,clear

n = 5;
h = 1/(n-2);
[xx,yy] = meshgrid( -h/2 : h : 1+h/2 );
x = xx(:);  y = yy(:);
[xV,yV] = meshgrid( 0:h:1, h/2:h:1-h/2 );
xV = xV(:);  yV = yV(:);
[xH,yH] = meshgrid( h/2:h:1-h/2, 0:h:1 );
xH = xH(:);  yH = yH(:);

idx = knnsearch( [x,y], [xV,yV], 'k', stencilSize );
