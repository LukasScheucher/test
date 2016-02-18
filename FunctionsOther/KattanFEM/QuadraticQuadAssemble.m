function y = QuadraticQuadAssemble(K,k,i,j,m,p,q,r,s,t)
%QuadraticQuadAssemble   This function assembles the element 
%                        stiffness matrix k of the quadratic 
%                        quadrilateral element with nodes i, j, 
%                        m, p, q, r, s, and t into the global 
%                        stiffness matrix K.
%                        This function returns the global 
%                        stiffness matrix K after the element 
%                        stiffness matrix k is assembled.
K(2*i-1,2*i-1) = K(2*i-1,2*i-1) + k(1,1);
K(2*i-1,2*i) = K(2*i-1,2*i) + k(1,2);
K(2*i-1,2*j-1) = K(2*i-1,2*j-1) + k(1,3);
K(2*i-1,2*j) = K(2*i-1,2*j) + k(1,4);
K(2*i-1,2*m-1) = K(2*i-1,2*m-1) + k(1,5);
K(2*i-1,2*m) = K(2*i-1,2*m) + k(1,6);
K(2*i-1,2*p-1) = K(2*i-1,2*p-1) + k(1,7);
K(2*i-1,2*p) = K(2*i-1,2*p) + k(1,8);
K(2*i-1,2*q-1) = K(2*i-1,2*q-1) + k(1,9);
K(2*i-1,2*q) = K(2*i-1,2*q) + k(1,10);
K(2*i-1,2*r-1) = K(2*i-1,2*r-1) + k(1,11);
K(2*i-1,2*r) = K(2*i-1,2*r) + k(1,12);
K(2*i-1,2*s-1) = K(2*i-1,2*s-1) + k(1,13);
K(2*i-1,2*s) = K(2*i-1,2*s) + k(1,14);
K(2*i-1,2*t-1) = K(2*i-1,2*t-1) + k(1,15);
K(2*i-1,2*t) = K(2*i-1,2*t) + k(1,16);
K(2*i,2*i-1) = K(2*i,2*i-1) + k(2,1);
K(2*i,2*i) = K(2*i,2*i) + k(2,2);
K(2*i,2*j-1) = K(2*i,2*j-1) + k(2,3);
K(2*i,2*j) = K(2*i,2*j) + k(2,4);
K(2*i,2*m-1) = K(2*i,2*m-1) + k(2,5);
K(2*i,2*m) = K(2*i,2*m) + k(2,6);
K(2*i,2*p-1) = K(2*i,2*p-1) + k(2,7);
K(2*i,2*p) = K(2*i,2*p) + k(2,8);
K(2*i,2*q-1) = K(2*i,2*q-1) + k(2,9);
K(2*i,2*q) = K(2*i,2*q) + k(2,10);
K(2*i,2*r-1) = K(2*i,2*r-1) + k(2,11);
K(2*i,2*r) = K(2*i,2*r) + k(2,12);
K(2*i,2*s-1) = K(2*i,2*s-1) + k(2,13);
K(2*i,2*s) = K(2*i,2*s) + k(2,14);
K(2*i,2*t-1) = K(2*i,2*t-1) + k(2,15);
K(2*i,2*t) = K(2*i,2*t) + k(2,16);
K(2*j-1,2*i-1) = K(2*j-1,2*i-1) + k(3,1);
K(2*j-1,2*i) = K(2*j-1,2*i) + k(3,2);
K(2*j-1,2*j-1) = K(2*j-1,2*j-1) + k(3,3);
K(2*j-1,2*j) = K(2*j-1,2*j) + k(3,4);
K(2*j-1,2*m-1) = K(2*j-1,2*m-1) + k(3,5);
K(2*j-1,2*m) = K(2*j-1,2*m) + k(3,6);
K(2*j-1,2*p-1) = K(2*j-1,2*p-1) + k(3,7);
K(2*j-1,2*p) = K(2*j-1,2*p) + k(3,8);
K(2*j-1,2*q-1) = K(2*j-1,2*q-1) + k(3,9);
K(2*j-1,2*q) = K(2*j-1,2*q) + k(3,10);
K(2*j-1,2*r-1) = K(2*j-1,2*r-1) + k(3,11);
K(2*j-1,2*r) = K(2*j-1,2*r) + k(3,12);
K(2*j-1,2*s-1) = K(2*j-1,2*s-1) + k(3,13);
K(2*j-1,2*s) = K(2*j-1,2*s) + k(3,14);
K(2*j-1,2*t-1) = K(2*j-1,2*t-1) + k(3,15);
K(2*j-1,2*t) = K(2*j-1,2*t) + k(3,16);
K(2*j,2*i-1) = K(2*j,2*i-1) + k(4,1);
K(2*j,2*i) = K(2*j,2*i) + k(4,2);
K(2*j,2*j-1) = K(2*j,2*j-1) + k(4,3);
K(2*j,2*j) = K(2*j,2*j) + k(4,4);
K(2*j,2*m-1) = K(2*j,2*m-1) + k(4,5);
K(2*j,2*m) = K(2*j,2*m) + k(4,6);
K(2*j,2*p-1) = K(2*j,2*p-1) + k(4,7);
K(2*j,2*p) = K(2*j,2*p) + k(4,8);
K(2*j,2*q-1) = K(2*j,2*q-1) + k(4,9);
K(2*j,2*q) = K(2*j,2*q) + k(4,10);
K(2*j,2*r-1) = K(2*j,2*r-1) + k(4,11);
K(2*j,2*r) = K(2*j,2*r) + k(4,12);
K(2*j,2*s-1) = K(2*j,2*s-1) + k(4,13);
K(2*j,2*s) = K(2*j,2*s) + k(4,14);
K(2*j,2*t-1) = K(2*j,2*t-1) + k(4,15);
K(2*j,2*t) = K(2*j,2*t) + k(4,16);
K(2*m-1,2*i-1) = K(2*m-1,2*i-1) + k(5,1);
K(2*m-1,2*i) = K(2*m-1,2*i) + k(5,2);
K(2*m-1,2*j-1) = K(2*m-1,2*j-1) + k(5,3);
K(2*m-1,2*j) = K(2*m-1,2*j) + k(5,4);
K(2*m-1,2*m-1) = K(2*m-1,2*m-1) + k(5,5);
K(2*m-1,2*m) = K(2*m-1,2*m) + k(5,6);
K(2*m-1,2*p-1) = K(2*m-1,2*p-1) + k(5,7);
K(2*m-1,2*p) = K(2*m-1,2*p) + k(5,8);
K(2*m-1,2*q-1) = K(2*m-1,2*q-1) + k(5,9);
K(2*m-1,2*q) = K(2*m-1,2*q) + k(5,10);
K(2*m-1,2*r-1) = K(2*m-1,2*r-1) + k(5,11);
K(2*m-1,2*r) = K(2*m-1,2*r) + k(5,12);
K(2*m-1,2*s-1) = K(2*m-1,2*s-1) + k(5,13);
K(2*m-1,2*s) = K(2*m-1,2*s) + k(5,14);
K(2*m-1,2*t-1) = K(2*m-1,2*t-1) + k(5,15);
K(2*m-1,2*t) = K(2*m-1,2*t) + k(5,16);
K(2*m,2*i-1) = K(2*m,2*i-1) + k(6,1);
K(2*m,2*i) = K(2*m,2*i) + k(6,2);
K(2*m,2*j-1) = K(2*m,2*j-1) + k(6,3);
K(2*m,2*j) = K(2*m,2*j) + k(6,4);
K(2*m,2*m-1) = K(2*m,2*m-1) + k(6,5);
K(2*m,2*m) = K(2*m,2*m) + k(6,6);
K(2*m,2*p-1) = K(2*m,2*p-1) + k(6,7);
K(2*m,2*p) = K(2*m,2*p) + k(6,8);
K(2*m,2*q-1) = K(2*m,2*q-1) + k(6,9);
K(2*m,2*q) = K(2*m,2*q) + k(6,10);
K(2*m,2*r-1) = K(2*m,2*r-1) + k(6,11);
K(2*m,2*r) = K(2*m,2*r) + k(6,12);
K(2*m,2*s-1) = K(2*m,2*s-1) + k(6,13);
K(2*m,2*s) = K(2*m,2*s) + k(6,14);
K(2*m,2*t-1) = K(2*m,2*t-1) + k(6,15);
K(2*m,2*t) = K(2*m,2*t) + k(6,16);
K(2*p-1,2*i-1) = K(2*p-1,2*i-1) + k(7,1);
K(2*p-1,2*i) = K(2*p-1,2*i) + k(7,2);
K(2*p-1,2*j-1) = K(2*p-1,2*j-1) + k(7,3);
K(2*p-1,2*j) = K(2*p-1,2*j) + k(7,4);
K(2*p-1,2*m-1) = K(2*p-1,2*m-1) + k(7,5);
K(2*p-1,2*m) = K(2*p-1,2*m) + k(7,6);
K(2*p-1,2*p-1) = K(2*p-1,2*p-1) + k(7,7);
K(2*p-1,2*p) = K(2*p-1,2*p) + k(7,8);
K(2*p-1,2*q-1) = K(2*p-1,2*q-1) + k(7,9);
K(2*p-1,2*q) = K(2*p-1,2*q) + k(7,10);
K(2*p-1,2*r-1) = K(2*p-1,2*r-1) + k(7,11);
K(2*p-1,2*r) = K(2*p-1,2*r) + k(7,12);
K(2*p-1,2*s-1) = K(2*p-1,2*s-1) + k(7,13);
K(2*p-1,2*s) = K(2*p-1,2*s) + k(7,14);
K(2*p-1,2*t-1) = K(2*p-1,2*t-1) + k(7,15);
K(2*p-1,2*t) = K(2*p-1,2*t) + k(7,16);
K(2*p,2*i-1) = K(2*p,2*i-1) + k(8,1);
K(2*p,2*i) = K(2*p,2*i) + k(8,2);
K(2*p,2*j-1) = K(2*p,2*j-1) + k(8,3);
K(2*p,2*j) = K(2*p,2*j) + k(8,4);
K(2*p,2*m-1) = K(2*p,2*m-1) + k(8,5);
K(2*p,2*m) = K(2*p,2*m) + k(8,6);
K(2*p,2*p-1) = K(2*p,2*p-1) + k(8,7);
K(2*p,2*p) = K(2*p,2*p) + k(8,8);
K(2*p,2*q-1) = K(2*p,2*q-1) + k(8,9);
K(2*p,2*q) = K(2*p,2*q) + k(8,10);
K(2*p,2*r-1) = K(2*p,2*r-1) + k(8,11);
K(2*p,2*r) = K(2*p,2*r) + k(8,12);
K(2*p,2*s-1) = K(2*p,2*s-1) + k(8,13);
K(2*p,2*s) = K(2*p,2*s) + k(8,14);
K(2*p,2*t-1) = K(2*p,2*t-1) + k(8,15);
K(2*p,2*t) = K(2*p,2*t) + k(8,16);
K(2*q-1,2*i-1) = K(2*q-1,2*i-1) + k(9,1);
K(2*q-1,2*i) = K(2*q-1,2*i) + k(9,2);
K(2*q-1,2*j-1) = K(2*q-1,2*j-1) + k(9,3);
K(2*q-1,2*j) = K(2*q-1,2*j) + k(9,4);
K(2*q-1,2*m-1) = K(2*q-1,2*m-1) + k(9,5);
K(2*q-1,2*m) = K(2*q-1,2*m) + k(9,6);
K(2*q-1,2*p-1) = K(2*q-1,2*p-1) + k(9,7);
K(2*q-1,2*p) = K(2*q-1,2*p) + k(9,8);
K(2*q-1,2*q-1) = K(2*q-1,2*q-1) + k(9,9);
K(2*q-1,2*q) = K(2*q-1,2*q) + k(9,10);
K(2*q-1,2*r-1) = K(2*q-1,2*r-1) + k(9,11);
K(2*q-1,2*r) = K(2*q-1,2*r) + k(9,12);
K(2*q-1,2*s-1) = K(2*q-1,2*s-1) + k(9,13);
K(2*q-1,2*s) = K(2*q-1,2*s) + k(9,14);
K(2*q-1,2*t-1) = K(2*q-1,2*t-1) + k(9,15);
K(2*q-1,2*t) = K(2*q-1,2*t) + k(9,16);
K(2*q,2*i-1) = K(2*q,2*i-1) + k(10,1);
K(2*q,2*i) = K(2*q,2*i) + k(10,2);
K(2*q,2*j-1) = K(2*q,2*j-1) + k(10,3);
K(2*q,2*j) = K(2*q,2*j) + k(10,4);
K(2*q,2*m-1) = K(2*q,2*m-1) + k(10,5);
K(2*q,2*m) = K(2*q,2*m) + k(10,6);
K(2*q,2*p-1) = K(2*q,2*p-1) + k(10,7);
K(2*q,2*p) = K(2*q,2*p) + k(10,8);
K(2*q,2*q-1) = K(2*q,2*q-1) + k(10,9);
K(2*q,2*q) = K(2*q,2*q) + k(10,10);
K(2*q,2*r-1) = K(2*q,2*r-1) + k(10,11);
K(2*q,2*r) = K(2*q,2*r) + k(10,12);
K(2*q,2*s-1) = K(2*q,2*s-1) + k(10,13);
K(2*q,2*s) = K(2*q,2*s) + k(10,14);
K(2*q,2*t-1) = K(2*q,2*t-1) + k(10,15);
K(2*q,2*t) = K(2*q,2*t) + k(10,16);
K(2*r-1,2*i-1) = K(2*r-1,2*i-1) + k(11,1);
K(2*r-1,2*i) = K(2*r-1,2*i) + k(11,2);
K(2*r-1,2*j-1) = K(2*r-1,2*j-1) + k(11,3);
K(2*r-1,2*j) = K(2*r-1,2*j) + k(11,4);
K(2*r-1,2*m-1) = K(2*r-1,2*m-1) + k(11,5);
K(2*r-1,2*m) = K(2*r-1,2*m) + k(11,6);
K(2*r-1,2*p-1) = K(2*r-1,2*p-1) + k(11,7);
K(2*r-1,2*p) = K(2*r-1,2*p) + k(11,8);
K(2*r-1,2*q-1) = K(2*r-1,2*q-1) + k(11,9);
K(2*r-1,2*q) = K(2*r-1,2*q) + k(11,10);
K(2*r-1,2*r-1) = K(2*r-1,2*r-1) + k(11,11);
K(2*r-1,2*r) = K(2*r-1,2*r) + k(11,12);
K(2*r-1,2*s-1) = K(2*r-1,2*s-1) + k(11,13);
K(2*r-1,2*s) = K(2*r-1,2*s) + k(11,14);
K(2*r-1,2*t-1) = K(2*r-1,2*t-1) + k(11,15);
K(2*r-1,2*t) = K(2*r-1,2*t) + k(11,16);
K(2*r,2*i-1) = K(2*r,2*i-1) + k(12,1);
K(2*r,2*i) = K(2*r,2*i) + k(12,2);
K(2*r,2*j-1) = K(2*r,2*j-1) + k(12,3);
K(2*r,2*j) = K(2*r,2*j) + k(12,4);
K(2*r,2*m-1) = K(2*r,2*m-1) + k(12,5);
K(2*r,2*m) = K(2*r,2*m) + k(12,6);
K(2*r,2*p-1) = K(2*r,2*p-1) + k(12,7);
K(2*r,2*p) = K(2*r,2*p) + k(12,8);
K(2*r,2*q-1) = K(2*r,2*q-1) + k(12,9);
K(2*r,2*q) = K(2*r,2*q) + k(12,10);
K(2*r,2*r-1) = K(2*r,2*r-1) + k(12,11);
K(2*r,2*r) = K(2*r,2*r) + k(12,12);
K(2*r,2*s-1) = K(2*r,2*s-1) + k(12,13);
K(2*r,2*s) = K(2*r,2*s) + k(12,14);
K(2*r,2*t-1) = K(2*r,2*t-1) + k(12,15);
K(2*r,2*t) = K(2*r,2*t) + k(12,16);
K(2*s-1,2*i-1) = K(2*s-1,2*i-1) + k(13,1);
K(2*s-1,2*i) = K(2*s-1,2*i) + k(13,2);
K(2*s-1,2*j-1) = K(2*s-1,2*j-1) + k(13,3);
K(2*s-1,2*j) = K(2*s-1,2*j) + k(13,4);
K(2*s-1,2*m-1) = K(2*s-1,2*m-1) + k(13,5);
K(2*s-1,2*m) = K(2*s-1,2*m) + k(13,6);
K(2*s-1,2*p-1) = K(2*s-1,2*p-1) + k(13,7);
K(2*s-1,2*p) = K(2*s-1,2*p) + k(13,8);
K(2*s-1,2*q-1) = K(2*s-1,2*q-1) + k(13,9);
K(2*s-1,2*q) = K(2*s-1,2*q) + k(13,10);
K(2*s-1,2*r-1) = K(2*s-1,2*r-1) + k(13,11);
K(2*s-1,2*r) = K(2*s-1,2*r) + k(13,12);
K(2*s-1,2*s-1) = K(2*s-1,2*s-1) + k(13,13);
K(2*s-1,2*s) = K(2*s-1,2*s) + k(13,14);
K(2*s-1,2*t-1) = K(2*s-1,2*t-1) + k(13,15);
K(2*s-1,2*t) = K(2*s-1,2*t) + k(13,16);
K(2*s,2*i-1) = K(2*s,2*i-1) + k(14,1);
K(2*s,2*i) = K(2*s,2*i) + k(14,2);
K(2*s,2*j-1) = K(2*s,2*j-1) + k(14,3);
K(2*s,2*j) = K(2*s,2*j) + k(14,4);
K(2*s,2*m-1) = K(2*s,2*m-1) + k(14,5);
K(2*s,2*m) = K(2*s,2*m) + k(14,6);
K(2*s,2*p-1) = K(2*s,2*p-1) + k(14,7);
K(2*s,2*p) = K(2*s,2*p) + k(14,8);
K(2*s,2*q-1) = K(2*s,2*q-1) + k(14,9);
K(2*s,2*q) = K(2*s,2*q) + k(14,10);
K(2*s,2*r-1) = K(2*s,2*r-1) + k(14,11);
K(2*s,2*r) = K(2*s,2*r) + k(14,12);
K(2*s,2*s-1) = K(2*s,2*s-1) + k(14,13);
K(2*s,2*s) = K(2*s,2*s) + k(14,14);
K(2*s,2*t-1) = K(2*s,2*t-1) + k(14,15);
K(2*s,2*t) = K(2*s,2*t) + k(14,16);
K(2*t-1,2*i-1) = K(2*t-1,2*i-1) + k(15,1);
K(2*t-1,2*i) = K(2*t-1,2*i) + k(15,2);
K(2*t-1,2*j-1) = K(2*t-1,2*j-1) + k(15,3);
K(2*t-1,2*j) = K(2*t-1,2*j) + k(15,4);
K(2*t-1,2*m-1) = K(2*t-1,2*m-1) + k(15,5);
K(2*t-1,2*m) = K(2*t-1,2*m) + k(15,6);
K(2*t-1,2*p-1) = K(2*t-1,2*p-1) + k(15,7);
K(2*t-1,2*p) = K(2*t-1,2*p) + k(15,8);
K(2*t-1,2*q-1) = K(2*t-1,2*q-1) + k(15,9);
K(2*t-1,2*q) = K(2*t-1,2*q) + k(15,10);
K(2*t-1,2*r-1) = K(2*t-1,2*r-1) + k(15,11);
K(2*t-1,2*r) = K(2*t-1,2*r) + k(15,12);
K(2*t-1,2*s-1) = K(2*t-1,2*s-1) + k(15,13);
K(2*t-1,2*s) = K(2*t-1,2*s) + k(15,14);
K(2*t-1,2*t-1) = K(2*t-1,2*t-1) + k(15,15);
K(2*t-1,2*t) = K(2*t-1,2*t) + k(15,16);
K(2*t,2*i-1) = K(2*t,2*i-1) + k(16,1);
K(2*t,2*i) = K(2*t,2*i) + k(16,2);
K(2*t,2*j-1) = K(2*t,2*j-1) + k(16,3);
K(2*t,2*j) = K(2*t,2*j) + k(16,4);
K(2*t,2*m-1) = K(2*t,2*m-1) + k(16,5);
K(2*t,2*m) = K(2*t,2*m) + k(16,6);
K(2*t,2*p-1) = K(2*t,2*p-1) + k(16,7);
K(2*t,2*p) = K(2*t,2*p) + k(16,8);
K(2*t,2*q-1) = K(2*t,2*q-1) + k(16,9);
K(2*t,2*q) = K(2*t,2*q) + k(16,10);
K(2*t,2*r-1) = K(2*t,2*r-1) + k(16,11);
K(2*t,2*r) = K(2*t,2*r) + k(16,12);
K(2*t,2*s-1) = K(2*t,2*s-1) + k(16,13);
K(2*t,2*s) = K(2*t,2*s) + k(16,14);
K(2*t,2*t-1) = K(2*t,2*t-1) + k(16,15);
K(2*t,2*t) = K(2*t,2*t) + k(16,16);
y = K;



