% Conservative transport scheme using PHS RBFs and polynomials for integration across cell walls.
% Greg Barnett
% August 9th, 2016

function [ psi, psiExact ] = fluxRBFFD( uvCase, h, k, rbforder, polyorder, stencilSize )

t = 0 : k : 5;                                          % vector of time values
rksteps = 3;
plotSolution = 1;

phi = @(x,y)  (x.^2+y.^2) .^ (rbforder/2);              % RBF which will be shifted to create RBF part of basis
n = 1/h;                                                % number of cells going across the domain
k = t(2) - t(1);                                        % time elapsed during one time step
[xx,yy] = meshgrid( -5*h/2 : h : 1+5*h/2 );

% psi = exp( -100*( (xx-.5).^2 + (yy-.5).^2 ) );
psi = exp( -625*( (xx-.5).^4 + (yy-.5).^4 ) );
psi = psi(:);											% initial condition for scalar field psi
psiExact = psi;

x = xx(:);  y = yy(:);                                  % location of cell-averaged values
% psiMaxBndry = psi( abs(x-.5)<h/4 & abs(y-1+h/2)<h/4 )
ii = x>0 & x<1 & y>0 & y<1;
[xV,yV] = meshgrid( 0:h:1, h/2:h:1-h/2 );
xV = xV(:);  yV = yV(:);                                % location of vertical cell walls (midpoint)
[xH,yH] = meshgrid( h/2:h:1-h/2, 0:h:1 );
xH = xH(:);  yH = yH(:);                                % location of horizontal cell walls (midpoint)

idxV = knnsearch( [x,y], [xV,yV], 'k', stencilSize, 'distance', 'chebychev' );
idxH = knnsearch( [x,y], [xH,yH], 'k', stencilSize, 'distance', 'chebychev' );

Wv = quadWeights( phi, @poly, h, x, y, idxV, xV, yV, 1, rbforder, polyorder );
Wh = quadWeights( phi, @poly, h, x, y, idxH, xH, yH, 0, rbforder, polyorder );

Wlr = Wv(1:n^2,:) - Wv(n+1:n^2+n,:);
specialIndex = [];
for i = 1 : n
	specialIndex = [ specialIndex, (i-1)*(n+1)+(1:n) ];
end
Wbt = Wh(specialIndex,:) - Wh(specialIndex+1,:);

[ xxtmp, yytmp ] = meshgrid( 0:h:1 );
for i = 1 : length(xV)
	figure(1),clf
		mesh( xxtmp, yytmp, -1*ones(size(xxtmp)) )
		view(2)
		hold('on')
		plot(x,y,'k.',xH(i),yH(i),'r.', x(idxH(i,:)), y(idxH(i,:)),'ro')
		axis('equal')
	drawnow,hold('off'),pause
end

U = @(t) u(x,y,t,uvCase);  V = @(t) v(x,y,t,uvCase);
f = @(t,psi)  odefun( t, psi, U(t), V(t), idxV, idxH, n, h, x, y, xV, yV, xH, yH, Wlr, Wbt, ii );

ep = 1e-10;
for i = 1 : length(t)-1
	psi = rk( t(i), psi, k, f, rksteps );
	if plotSolution == 1
		if abs( round(t(i+1)*100) - t(i+1)*100 ) <= ep && mod( round(t(i+1)*100), 5) == 0
			psi = reshape( psi, sqrt(length(psi)), sqrt(length(psi)) );
			psi([1,2,3,end,end-1,end-2],:) = psi([end-5,end-4,end-3,6,5,4],:);
			psi(:,[1,2,3,end,end-1,end-2]) = psi(:,[end-5,end-4,end-3,6,5,4]);
			figure(1),clf
				[~,H]=contourf( xx, yy, psi, -.05:.1:1.05 );  set( H, 'lineStyle', 'none' )
				% surf( xx, yy, psi ),view(2),shading('interp'),lighting('phong')
				% contour( xx, yy, psi, -.05:.1:.95 )
				axis( 'equal', [0,1,0,1] )
				set( gca, 'xTickLabel', [], 'yTickLabel', [] )
				caxis( [-.05,1.05] )
				colorbar
				colormap( parula(11) )
				title(sprintf('t=%g, min=%g, max=%g, mass=%g',t(i+1),min(min(psi)),max(max(psi)),h^2*sum(psi(ii))))
				% title(sprintf('RBF=r^%g, poly=%g, stencil=%g',rbforder,polyorder,stencilSize))
			drawnow
			psi = psi(:);
		end
	end
end

psi = reshape( psi, sqrt(length(psi)), sqrt(length(psi)) );
psi([1,2,3,end,end-1,end-2],:) = psi([end-5,end-4,end-3,6,5,4],:);
psi(:,[1,2,3,end,end-1,end-2]) = psi(:,[end-5,end-4,end-3,6,5,4]);
psi = psi(:);
h = 1/101;
[X,Y] = meshgrid( h/2 : h : 1-h/2 );
X = X(:);  Y = Y(:);
idx = knnsearch( [x,y], [X,Y], 'k', 1 );
idx = knnsearch( [x,y], [x(idx),y(idx)], 'k', 25 );
W = rbffd( idx, x, y, X, Y, @(x,y)(x.^2+y.^2).^(5/2), @(x,y)(x.^2+y.^2).^(5/2), [1], 3 );
psi = W * psi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THE ODE FUNCTION
% U and V should be functions of t only.

function z = odefun( t, psi, U, V, idxV, idxH, n, h, x, y, xV, yV, xH, yH, Wlr, Wbt, ii )

psi = reshape( psi, n+6, n+6 );
psi( [1,2,3,end,end-1,end-2], : ) = psi( [end-5,end-4,end-3,6,5,4], : );
psi( :, [1,2,3,end,end-1,end-2] ) = psi( :, [end-5,end-4,end-3,6,5,4] );
psi = psi(:);

z = zeros( size(psi) );
z(ii) = 1/h^2 * ( Wlr*(psi.*U) + Wbt*(psi.*V) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUADRATURE WEIGHTS FUNCTION

function W = quadWeights( phi, poly, h, x, y, idx, xe, ye, vert, rbforder, polyorder )
eps = 1e-12;
stencilSize = size(idx,2);  ne = length(xe);
X = x(idx);  Y = y(idx);
if rbforder == 1
	antiderivative = @(th) tan(th);
elseif rbforder == 3
	antiderivative = @(th) 1/3 * ( cos(2*th) + 2 ) * tan(th) * sec(th)^2;
elseif rbforder == 5
	antiderivative = @(th) 1/15 * ( 6*cos(2*th) + cos(4*th) + 8 ) * tan(th) * sec(th)^4;
elseif rbforder == 7
	antiderivative = @(th) 1/70 * ( 29*cos(2*th) + 8*cos(4*th) + cos(6*th) + 32 ) * tan(th) * sec(th)^6;
elseif rbforder == 9
	antiderivative = @(th) 1/315 * ( 130*cos(2*th) + 46*cos(4*th) + 10*cos(6*th) + cos(8*th) + 128 ) * tan(th) * sec(th)^8;
end
W = zeros( ne, stencilSize );
for i = 1 : ne
    xn = X(i,:) - xe(i);  yn = Y(i,:) - ye(i);
    xx = meshgrid(xn);  yy = meshgrid(yn);
    P = poly(xn,yn,polyorder);
    A = [ phi(xx.'-xx,yy.'-yy), P.';  P, zeros(size(P,1),size(P,1)) ];
	b = zeros( 1, stencilSize );
    if vert == 1
		for j = 1 : length(xn)
			if ( xn(j)>0 && yn(j)>0 ) || ( xn(j)<0 && yn(j)<0 )
				th1 = atan( ( yn(j) - h/2 ) / xn(j) );
				th2 = atan( ( yn(j) + h/2 ) / xn(j) );
			elseif ( xn(j)<0 && yn(j)>0 ) || ( xn(j)>0 && yn(j) < 0 )
				th1 = atan( ( yn(j) - h/2 ) / -xn(j) );
				th2 = atan( ( yn(j) + h/2 ) / -xn(j) );
			else
				th1 = atan( -h/2 / xn(j) );
				th2 = atan( h/2 / xn(j) );
			end
			b(j) = xn(j)^(rbforder+1) * abs( antiderivative(th2) - antiderivative(th1) );
		end
		if polyorder == 1
			b = [ b, [h,0,0] ];
		elseif polyorder == 2
			b = [ b, [h,0,0], [0,0,h^3/12] ];
		elseif polyorder == 3
			b = [ b, [h,0,0], [0,0,h^3/12], [0,0,0,0] ];
		elseif polyorder == 4
			b = [ b, [h,0,0], [0,0,h^3/12], [0,0,0,0], [0,0,0,0,h^5/80] ];
		end
    else
		for j = 1 : length(xn)
			if ( xn(j)>0 && yn(j)>0 ) || ( xn(j)<0 && yn(j)<0 )
				th1 = atan( ( xn(j) + h/2 ) / yn(j) );
				th2 = atan( ( xn(j) - h/2 ) / yn(j) );
			elseif ( xn(j)<0 && yn(j)>0 ) || ( xn(j)>0 && yn(j) < 0 )
				th1 = atan( ( xn(j) - h/2 ) / -yn(j) );
				th2 = atan( ( xn(j) + h/2 ) / -yn(j) );
			else
				th1 = atan( -h/2 / yn(j) );
				th2 = atan( h/2 / yn(j) );
			end
			b(j) = yn(j)^(rbforder+1) * abs( antiderivative(th2) - antiderivative(th1) );
		end
		if polyorder == 1
			b = [ b, [h,0,0] ];
		elseif polyorder == 2
			b = [ b, [h,0,0], [h^3/12,0,0] ];
		elseif polyorder == 3
			b = [ b, [h,0,0], [h^3/12,0,0], [0,0,0,0] ];
		elseif polyorder == 4
			b = [ b, [h,0,0], [h^3/12,0,0], [0,0,0,0], [h^5/80,0,0,0,0] ];
		end
    end
	invA = inv(A),pause
    w = b / A;
    W(i,:) = w( 1 : stencilSize );
end
ii = repmat( 1:ne, stencilSize, 1 );
W = sparse( ii, idx.', W.', ne, length(x), ne*stencilSize );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%HORIZONTAL VELOCITY

function z = u(x,y,t,uvCase)
if uvCase == 1
	z = ones(size(x));
elseif uvCase == 2
	z = ( 1 + cos(2*pi*x).*sin(2*pi*y) ) ./ 2;
elseif uvCase == 3
	z = cos(2*pi*y);
elseif uvCase == 4
	z = ones(size(x));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VERTICAL VELOCITY

function z = v(x,y,t,uvCase)
if uvCase == 1
	z = cos(2*pi*x) .* sin(2*pi*y);
elseif uvCase == 2
	z = ones(size(x));
elseif uvCase == 3
	z = ones(size(x));
elseif uvCase == 4
	z = ones(size(x));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%RUNGE KUTTA

function U = rk(t, U, k, odeFun, rkSteps)
if rkSteps == 1
    U = U + k*odeFun(t,U);
elseif rkSteps == 2
    s1 = odeFun(t,U);                                 %2nd order Runge-Kutta, step1
    s1 = odeFun(t+k/2,U+k/2*s1);                      %step 2
    U  = U + k*s1;                                    %Combine for new value
elseif rkSteps == 3
    s1 = odeFun(t,U);                                 %3rd order Runge-Kutta, step 1
    s2 = odeFun(t+k/3,U+k/3*s1);                      %step 2
    s2 = odeFun(t+2*k/3,U+2*k/3*s2);                  %step 3
    U  = U + k/4 * ( s1 + 3*s2 );                     %Combine for new value
else
    s1 = odeFun( t, U );                              %4th order Runge-Kutta, step 1
    s2 = odeFun( t+k/2, U+k/2*s1 );                   %step 2
    s3 = odeFun( t+k/2, U+k/2*s2 );                   %step 3
    s4 = odeFun( t+k, U+k*s3 );                       %step 4
    U  = U + k/6 * ( s1 + 2*s2 + 2*s3 + s4 );         %Combine for new value
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% JUST FOR INTERPOLATION TO FINER GRID FOR COMPARISON WITH PS RESULT

function w = rbffd( idx, x, z, xe, ze, phi, Lphi, v, polyorder )

ne = length( xe );
n = size( idx, 2 );
X = x(idx);  Z = z(idx);
w = zeros( n, ne );
ii = repmat( 1:ne, n, 1 );

parfor i = 1:ne
    w(:,i) = polyharmonicWeights2( X(i,:), Z(i,:), xe(i), ze(i), ...
        phi, Lphi, v, polyorder );
end

w = sparse( ii, idx.', w, ne, length(x), n*ne );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% JUST FOR INTERPOLATION TO FINER GRID FOR COMPARISON WITH PS RESULT

function w = polyharmonicWeights2( x, z, xe, ze, phi, Lphi, v, polyorder )

x = x - xe;  z = z - ze;
p = poly( x, z, polyorder );
B = [ Lphi(-x,-z), v, zeros(1,size(p,1)-length(v)) ];
x = meshgrid(x);  z = meshgrid(z);
A = [ phi(x.'-x,z.'-z), p.'; p, zeros( size(p,1), size(p,1) ) ];
w = B / A;
w = w( 1:length(x) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=poly(x,z,order)

if order==-1
    y=[];
elseif order==0
    y=ones(1,length(x));
elseif order==1
    y=[ ...
     ones(1,length(x)); ...
     x;z
     ];
elseif order==2
    y=[ ...
     ones(1,length(x)); ...
     x;z; ...
     x.^2;x.*z;z.^2
     ];
elseif order==3
    y=[ ...
     ones(1,length(x)); ...
     x;z; ...
     x.^2;x.*z;z.^2; ...
     x.^3;x.^2.*z;x.*z.^2;z.^3
     ];
elseif order==4
    y=[ ...
     ones(1,length(x)); ...
     x;z; ...
     x.^2;x.*z;z.^2; ...
     x.^3;x.^2.*z;x.*z.^2;z.^3; ...
     x.^4;x.^3.*z;x.^2.*z.^2;x.*z.^3;z.^4
     ];
elseif order==5
    y=[ ...
     ones(1,length(x)); ...
     x;z; ...
     x.^2;x.*z;z.^2; ...
     x.^3;x.^2.*z;x.*z.^2;z.^3; ...
     x.^4;x.^3.*z;x.^2.*z.^2;x.*z.^3;z.^4; ...
     x.^5;x.^4.*z;x.^3.*z.^2;x.^2.*z.^3;x.*z.^4;z.^5
     ];
elseif order==6
    y=[ ...
     ones(1,length(x)); ...
     x;z; ...
     x.^2;x.*z;z.^2; ...
     x.^3;x.^2.*z;x.*z.^2;z.^3; ...
     x.^4;x.^3.*z;x.^2.*z.^2;x.*z.^3;z.^4; ...
     x.^5;x.^4.*z;x.^3.*z.^2;x.^2.*z.^3;x.*z.^4;z.^5; ...
     x.^6;x.^5.*z;x.^4.*z.^2;x.^3.*z.^3;x.^2.*z.^4;x.*z.^5;z.^6
     ];
elseif order==7
     y=[ ...
     ones(1,length(x)); ...
     x;z; ...
     x.^2;x.*z;z.^2; ...
     x.^3;x.^2.*z;x.*z.^2;z.^3; ...
     x.^4;x.^3.*z;x.^2.*z.^2;x.*z.^3;z.^4; ...
     x.^5;x.^4.*z;x.^3.*z.^2;x.^2.*z.^3;x.*z.^4;z.^5; ...
     x.^6;x.^5.*z;x.^4.*z.^2;x.^3.*z.^3;x.^2.*z.^4;x.*z.^5;z.^6; ...
     x.^7;x.^6.*z;x.^5.*z.^2;x.^4.*z.^3;x.^3.*z.^4;x.^2.*z.^5;x.*z.^6;z.^7
     ];
elseif order==8
     y=[ ...
     ones(1,length(x)); ...
     x;z; ...
     x.^2;x.*z;z.^2; ...
     x.^3;x.^2.*z;x.*z.^2;z.^3; ...
     x.^4;x.^3.*z;x.^2.*z.^2;x.*z.^3;z.^4; ...
     x.^5;x.^4.*z;x.^3.*z.^2;x.^2.*z.^3;x.*z.^4;z.^5; ...
     x.^6;x.^5.*z;x.^4.*z.^2;x.^3.*z.^3;x.^2.*z.^4;x.*z.^5;z.^6; ...
     x.^7;x.^6.*z;x.^5.*z.^2;x.^4.*z.^3;x.^3.*z.^4;x.^2.*z.^5;x.*z.^6;z.^7; ...
     x.^8;x.^7.*z;x.^6.*z.^2;x.^5.*z.^3;x.^4.*z.^4;x.^3.*z.^5;x.^2.*z.^6;x.*z.^7;z.^8
     ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




























