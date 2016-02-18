function w = LinearBrickElementVolume(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8)
%LinearBrickElementVolume   This function returns the volume
%                           of the linear brick element
%                           whose first node has coordinates
%                           (x1,y1,z1), second node has 
%                           coordinates (x2,y2,z2), third node 
%                           has coordinates (x3,y3,z3), 
%                           fourth node has coordiantes
%                           (x4,y4,z4), fifth node has coordiantes
%                           (x5,y5,z5), sixth node has coordiantes
%                           (x6,y6,z6), seventh node has coordiantes
%                           (x7,y7,z7), and eighth node has 
%                           coordiantes (x8,y8,z8).
syms s t u;
N1 = (1-s)*(1-t)*(1+u)/8;
N2 = (1-s)*(1-t)*(1-u)/8;
N3 = (1-s)*(1+t)*(1-u)/8;
N4 = (1-s)*(1+t)*(1+u)/8;
N5 = (1+s)*(1-t)*(1+u)/8;
N6 = (1+s)*(1-t)*(1-u)/8;
N7 = (1+s)*(1+t)*(1-u)/8;
N8 = (1+s)*(1+t)*(1+u)/8;
x = N1*x1 + N2*x2 + N3*x3 + N4*x4 + N5*x5 + N6*x6 + N7*x7 + N8*x8;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4 + N5*y5 + N6*y6 + N7*y7 + N8*y8;
z = N1*z1 + N2*z2 + N3*z3 + N4*z4 + N5*z5 + N6*z6 + N7*z7 + N8*z8;
xs = diff(x,s);
xt = diff(x,t);
xu = diff(x,u);
ys = diff(y,s);
yt = diff(y,t);
yu = diff(y,u);
zs = diff(z,s);
zt = diff(z,t);
zu = diff(z,u);
J = xs*(yt*zu - zt*yu) - ys*(xt*zu - zt*xu) + zs*(xt*yu - yt*xu);
Jnew = simplify(J);
r = int(int(int(Jnew, u, -1, 1), t, -1, 1), s, -1, 1);
w = double(r);
