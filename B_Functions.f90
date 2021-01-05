Function Pressure(X,Y)
  Pressure = x**2+y**2
End Function

function RLinearInterp(d1,d2,x1,x2)
   RLinearInterp=(x1*d2+x2*d1)/(d1+d2)
end function
   
subroutine CalcGradPExact(X,Y,GPE)
REAL X,Y,GPE(2)
   GPE(1) = 2.0*x
   GPE(2) = 2.0*y
end subroutine

subroutine Velocity(x,y,V)
real x,y,V(2)
     V(1) = (1.0+x)*(1.0+y)!1.0+X
	 V(2) = x*y!1.0+Y
end subroutine

function DivVelocityExact(x,y,mode)
   real x,y
   integer mode
   if (mode.EQ.0) then 
                  DivVelocityExact = 1+x+y !2.0
   else 
                  DivVelocityExact = X**2+4*X*Y+2*X+Y**2+2*Y+1
   endif
end function

function RotVelocityExact(x,y)
   real x,y
   RotVelocityExact=1+x-y
end function

function LaplacianPExact(x,y)
   real x,y
   LaplacianPExact=4.0
end function
