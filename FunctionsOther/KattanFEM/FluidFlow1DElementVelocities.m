function y = FluidFlow1DElementVelocities(Kxx,L,p)
%FluidFlow1DElementVelocities   This function returns the element 
%                               velocity given the element   
%                               permeability coefficient Kxx, the 
%                               element length L, and the element 
%                               nodal potential (fluid head) 
%                               vector p.
y = -Kxx * [-1/L 1/L] * p;

