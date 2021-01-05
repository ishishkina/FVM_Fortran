Subroutine B_OutputFields(IO,NI,NJ,X,Y,P,GradP,GradPError,V,divV,divVError,rotV,rotVError, LapP,LapPError)
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::P
  Real,Dimension(0:NI,0:NJ,2)::GradP
  Real,Dimension(0:NI,0:NJ,2)::GradPError
  Real,Dimension(0:NI,0:NJ,2)::V
  Real,Dimension(0:NI,0:NJ)::divV
  Real,Dimension(0:NI,0:NJ)::divVError
  Real,Dimension(0:NI,0:NJ)::rotV
  Real,Dimension(0:NI,0:NJ)::rotVError
  Real,Dimension(0:NI,0:NJ)::LapP
  Real,Dimension(0:NI,0:NJ)::LapPError
  

  Write(IO,*) 'VARIABLES = "X", "Y", "P", "GradPx","GradPy","GradPErrorX","GradPErrorY", "Vx", "Vy"'
    Write(IO,*) ' ,  "divV", "divVError","rotV", "rotVError", "LapP", "LapPError" '
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-15]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ) 
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)
  Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)
  write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,1)
  write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,2)
  write(IO,'(100F14.7)') GradPError(1:NI-1,1:NJ-1,1)
  write(IO,'(100F14.7)') GradPError(1:NI-1,1:NJ-1,2)
  write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,1)
  write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,2)
  write(IO,'(100F14.7)') divV(1:NI-1,1:NJ-1)
  write(IO,'(100F14.7)') divVError(1:NI-1,1:NJ-1)
  write(IO,'(100F14.7)') rotV(1:NI-1,1:NJ-1)
  write(IO,'(100F14.7)') rotVError(1:NI-1,1:NJ-1)
  write(IO,'(100F14.7)') LapP(1:NI-1,1:NJ-1)
  write(IO,'(100F14.7)') LapPError(1:NI-1,1:NJ-1)
  
End Subroutine 
