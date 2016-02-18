function y = FluidFlow1DElementVFR(Kxx,L,p,A)
%FluidFlow1DElementVFR   This function returns the 
%                        element volumetric flow rate given the  
%                        element permeability coefficient Kxx, the 
%                        element length L, the element 
%                        nodal potential (fluid head) vector
%                        p, and the element cross-sectional
%                        area A.
y = -Kxx * [-1/L 1/L] * p *A;

