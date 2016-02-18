function y = FluidFlow1DElementStiffness(Kxx,A,L)
%FluidFlow1DElementStiffness   This function returns the element 
%                              stiffness matrix for a fluid flow 
%                              1D element with coefficient of 
%                              permeability Kxx, cross-sectional 
%                              area A, and length L. The size of 
%                              the element stiffness matrix is 
%                              2 x 2.
y = [Kxx*A/L -Kxx*A/L ; -Kxx*A/L Kxx*A/L];


