Subroutine B_CalcDivergence(mode,NI,NJ,V,divV,P,GradP,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)

INTEGER I,J, NI,NJ,IN,JN,mode
REAL CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),RN(2)
REAL IFaceCenter(NI,NJ-1,2),IFaceVector(NI,NJ-1,2),JFaceCenter(NI-1,NJ,2),JFaceVector(NI-1,NJ,2)
REAL divV(0:NI,0:NJ), V(0:NI,0:NJ,2), NCELL(4,2),RF(4,2),SF(4,2),VOL,RC(2),DC,DN,VF(2),PF, P(0:NI,0:NJ), GradP(0:NI,0:NJ,2)
     
     
     DO I=1,NI-1
       DO J=1, NJ-1
          
          NCELL(1,:)=[I-1,J]
          NCELL(2,:)=[I+1,J]
          NCELL(3,:)=[I,J-1]
          NCELL(4,:)=[I,J+1]	

          RF(1,:)=IFaceCenter(I,J,:)
		  RF(2,:)=IFaceCenter(I+1,J,:)
		  RF(3,:)=JFaceCenter(I,J,:)
		  RF(4,:)=JFaceCenter(I,J+1,:)
		  
		  SF(1,:)=-IFaceVector(I,J,:)
		  SF(2,:)=IFaceVector(I+1,J,:)
		  SF(3,:)=-JFaceVector(I,J,:)
		  SF(4,:)=JFaceVector(I,J+1,:)
		  
		  VOL=CellVolume(I,J)
		  RC(:)=CellCenter(I,J,:)
		  
		  DO IFACE=1,4
		     IN=NCELL(IFACE,1)
			 JN=NCELL(IFACE,2)
			 RN(:)=CellCenter(IN,JN,:)
			 DC=Norm2(RF(IFace,:)-RC(:))
			 DN=Norm2(RF(IFace,:)-RN(:))
			 
			 VF(1)=RLinearInterp(DC,DN,V(I,J,1),V(IN,JN,1))
			 VF(2)=RLinearInterp(DC,DN,V(I,J,2),V(IN,JN,2))
			 
			 SELECT CASE(mode)
			 CASE(0)
			      divV(I,J) = divV(I,J)+DOT_PRODUCT(VF(:),SF(IFACE,:))
				  
			 CASE(1)
			      PF = RLinearInterp(DC,DN,P(I,J),P(IN,JN))
			      divV(I,J) = divV(I,J)+DOT_PRODUCT(PF*VF(:),SF(IFACE,:))
				  
			 CASE(2)
			      IF (DOT_PRODUCT(SF(IFACE,:),VF(:)).GE.0.0)THEN
				      PF=P(I,J)
				  ELSE 
				      PF=P(IN,JN)
					  IF (DN.LT.1E-6) PF=2*P(IN,JN)-P(I,J)
				  END IF
				  divV(I,J) = divV(I,J)+DOT_PRODUCT(PF*VF(:),SF(IFACE,:))
				  
			 CASE(3)
			      IF (DOT_PRODUCT(SF(IFACE,:),VF(:)).GE.0.0)THEN
				      PF=P(I,J)+DOT_PRODUCT(RF(IFACE,:)-CellCenter(I,J,:),GradP(I,J,:))
				  ELSE 
				      PF=P(IN,JN)+DOT_PRODUCT(RF(IFACE,:)-CellCenter(IN,JN,:),GradP(IN,JN,:))
					  IF (DN.LT.1E-6) THEN
					      PN=2*P(IN,JN)-P(I,J)
						  GC=DOT_PRODUCT(GradP(I,J,:),CellCenter(I,J,:)-RF(IFACE,:))
						  GB=P(I,J)-P(IN,JN)
						  GN=4*GB-3*GC
						  PF=PN+GN
					  END IF
				  END IF
				  divV(I,J) = divV(I,J)+DOT_PRODUCT(PF*VF(:),SF(IFACE,:))
				  
			 END SELECT 
			 
			 

		  END DO
		  divV(I,J)=divV(I,J)/VOL
        END DO
     END DO
     
End Subroutine 
