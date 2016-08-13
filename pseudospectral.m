function [ psi, psiExact ] = pseudospectral( uvCase, h, k )			% use h = 1 / (odd#)

t = 0 : k : 5;
rksteps = 4;
plotSolution = 1;

x = h/2 : h : 1-h/2;
y = x;
[xx,yy] = meshgrid( x, y );
n = 1/h;
% psi = exp( -100*( (xx-.5).^2 + (yy-.5).^2 ) );
psi = exp( -625*( (xx-.5).^4 + (yy-.5).^4 ) );
psiExact = psi;
% psiMaxBndry = psi( abs(xx-.5)<h/4 & abs(yy-1+h/2)<h/4 )

Acos = zeros( n, (n-1)/2+1 );
Bcos = zeros( n, (n-1)/2+1 );
for j = 1 : (n-1)/2+1
	Acos(:,j) = cos( (j-1)*2*pi*y );
	Bcos(:,j) = -(j-1)*2*pi * sin( (j-1)*2*pi*y );
end
Asin = zeros( n, (n-1)/2 );
Bsin = zeros( n, (n-1)/2 );
for j = 1 : (n-1)/2
	Asin(:,j) = sin( j*2*pi*y );
	Bsin(:,j) = j*2*pi * cos( j*2*pi*y );
end
A = [ Acos, Asin ];
B = [ Bcos, Bsin ];
Wy = B / A;
Wx = Wy.';

U = @(t) u(xx,yy,t,uvCase);  V = @(t) v(xx,yy,t,uvCase);
f = @(t,psi)  odefun( t, psi, U(t), V(t), Wx, Wy );

ep = 1e-10;
for i = 1 : length(t)-1
	psi = rk( t(i), psi, k, f, rksteps );
	if plotSolution == 1
		if abs( round(t(i+1)*100) - t(i+1)*100 ) <= ep && mod( round(t(i+1)*100), 5) == 0
			figure(1),clf
				[~,H]=contourf( xx, yy, psi, -.05:.1:.95 );  set( H, 'lineStyle', 'none' );
				% surf( xx, yy, psi ),view(2),shading('interp'),lighting('phong')
				% contour( xx, yy, psi, -.05:.1:.95 )
				axis( 'equal', [0,1,0,1] )
				set( gca, 'xTickLabel', [], 'yTickLabel', [] )
				caxis( [-.05,1.05] )
				colorbar
				colormap( parula(11) )
				% set(gca,'fontSize',20,'xTick',[],'yTick',[]),title(sprintf('t = %g',t(i+1)),'fontSize',25)
				title(sprintf('t=%g, min=%g, max=%g, mass=%g',t(i+1),min(min(psi)),max(max(psi)),h^2*sum(sum(psi))))
			drawnow
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = odefun( t, psi, U, V, Wx, Wy )
z = -(psi.*U) * Wx - Wy * (psi.*V);

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
