Program Main

  character(*), parameter:: InputFile='input.txt',OutputFile='data.plt' ! names of input and output files
  character MeshFile*30, SolutionFile*30, ctmp      ! name of file with computational mesh
  integer, parameter:: IO = 12 ! input-output unit
  integer Niter, mode
  real,allocatable,dimension(:,:):: X,Y,P,CellVolume, divV, divVExact,divVError, rotV, rotVExact, rotVError ! scalar arrays
  real,allocatable,dimension(:,:,:)::GradP, GradPExact, GradPError, V
  real,allocatable,dimension(:,:)::LapP, LapPExact, LapPError
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector ! vector arrays

!===  READ INPUT FILE ===
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
  READ(IO,*) MeshFile  ! read name of file with computational mesh
  READ(IO,*) SolutionFile  ! read name of file with solution data
  CLOSE(IO)
  mode=0
!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile
  OPEN(IO,FILE = MeshFile)
  READ(IO,*) NI,NJ
  WRITE(*,*) 'NI, NJ = ',NI,NJ

!=== ALLOCATE ALL ARRAYS ===
  WRITE(*,*) 'Allocate arrays'       
  allocate(X(NI,NJ)) ! mesh nodes X-coordinates
  allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(GradP(0:NI,0:NJ,2))
  allocate(GradPExact(0:NI,0:NJ,2))
  allocate(GradPError(0:NI,0:NJ,2))
  allocate(V(0:NI,0:NJ,2))
  allocate(divV(0:NI,0:NJ))
  allocate(divVExact(0:NI,0:NJ))
  allocate(divVError(0:NI,0:NJ))
  allocate(rotV(0:NI,0:NJ))
  allocate(rotVExact(0:NI,0:NJ))
  allocate(rotVError(0:NI,0:NJ))
  allocate(LapP(0:NI,0:NJ))
  allocate(LapPExact(0:NI,0:NJ))
  allocate(LapPError(0:NI,0:NJ))
  allocate(CellVolume(NI-1,NJ-1))   ! Cell Volumes    
  allocate(CellCenter(0:NI,0:NJ,2)) ! Cell Centers
  allocate(IFaceCenter( NI,NJ-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector( NI,NJ-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter( NI-1,NJ,2)) ! Face Centers for J-faces
  allocate(JFaceVector( NI-1,NJ,2)) ! Face Vectors for I-faces

!===  READ GRID ===
  WRITE(*,*) 'Read mesh from file: ', MeshFile
  READ(IO,*) ((X(I,J),Y(I,J),rtmp,I=1,NI),J=1,NJ)
  !READ(IO,*) Iter
  CLOSE(IO)
  Niter=15
!=== CALCULATE METRIC ===
  WRITE(*,*) 'Calculate metric'       
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector) 
  
!=== READ SOLUTION ===
   WRITE(*,*) 'Read solution data from file: ', SolutionFile
   OPEN(IO,FILE=SolutionFile) 
   READ(IO,*) ctmp
   READ(IO,*) ctmp
   READ(IO,*) ((rtmp,rtmp,V(I,J,1),V(I,J,2),rtmp,P(I,J),rtmp,rtmp,I=0,NI),J=0,NJ)
   
!=== INITIATE FIELDS ===
  WRITE(*,*) 'Initiate fields'       
  DO  J = 0,NJ
    DO  I = 0,NI
!      P(I,J) = Pressure(CellCenter(I,J,1),CellCenter(I,J,2))
!	  CALL CalcGradPExact(CellCenter(I,J,1),CellCenter(I,J,2),GradPExact(I,J,:))
!	  CALL Velocity(CellCenter(I,J,1),CellCenter(I,J,2),V(I,J,:))
!	  divVExact(I,J) = DivVelocityExact(CellCenter(I,J,1),CellCenter(I,J,2),mode)
!	  rotVExact(I,J) = RotVelocityExact(CellCenter(I,J,1),CellCenter(I,J,2))
!	  LapPExact(I,J) = LaplacianPExact(CellCenter(I,J,1),CellCenter(I,J,2))
    ENDDO
  ENDDO

!=== CALCULATE GRADIENT ===
  GradP=0.0
  WRITE(*,*) 'Calculate gradient' 
  DO I=1,Niter  
  Call B_CalcGradient(NI,NJ,P,GradP,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)
!  GradPError=ABS(GradPExact-GradP)/GradPExact
 ! write(*,*) 'Maximum GradP-error:', maxval(GradPError(1:NI-1,1:NJ-1,:))
  END DO
  
 ! write(*,*) 'Maximum GradPx-error:', maxval(GradPError(1:NI-1,1:NJ-1,1))
 ! write(*,*) 'Maximum GradPy-error:', maxval(GradPError(1:NI-1,1:NJ-1,2))
!=== CALCULATE Divergence ===
  divV=0.0
  WRITE(*,*) 'Calculate divergence'  
  Call B_CalcDivergence(mode, NI,NJ,V,divV,P,GradP,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)
!  divVError=ABS(divVExact-divV)/divVExact
!  write(*,*) 'Maximum divV-error:', maxval(divVError(1:NI-1,1:NJ-1))
  
!=== CALCULATE Rotor ===
  rotV=0.0
  WRITE(*,*) 'Calculate rotor'  
  Call B_CalcRotor(NI,NJ,V,rotV, CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)
!  rotVError=ABS(rotVExact-rotV)/rotVExact
!  write(*,*) 'Maximum rotV-error:', maxval(rotVError(1:NI-1,1:NJ-1))
  
  !=== CALCULATE Laplacian ===
  LapP=0.0
  WRITE(*,*) 'Calculate laplacian'  
  Call B_CalcLaplacian(NI,NJ,P,GradP,LapP, CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)
 ! LapPError=ABS(LapPExact-LapP)/LapPExact
 ! write(*,*) 'Maximum LapP-error:', maxval(LapPError(1:NI-1,1:NJ-1))
  
!=== OUTPUT FIELDS ===
  WRITE(*,*) 'Output fields to file: ', OutputFile       
  Open(IO,FILE=OutputFile)
  Call B_OutputFields(IO,NI,NJ,X,Y,P,GradP,GradPError,V,divV,divVError,rotV,rotVError, LapP, LapPError)
  Close(IO)

END PROGRAM Main  
