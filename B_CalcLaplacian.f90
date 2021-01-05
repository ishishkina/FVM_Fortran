Subroutine B_CalcLaplacian(NI,NJ,P,GradP,LapP, CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)

INTEGER I,J, NI,NJ,IN,JN,mode
REAL CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),RN(2)
REAL IFaceCenter(NI,NJ-1,2),IFaceVector(NI,NJ-1,2),JFaceCenter(NI-1,NJ,2),JFaceVector(NI-1,NJ,2)
REAL NCELL(4,2),RF(4,2),SF(4,2),VOL,RC(2),DC,DN, DNC, VF(2),PF, P(0:NI,0:NJ), GradP(0:NI,0:NJ,2),LapP(0:NI,0:NJ), NF(2),dpdn,dodn_c
REAL RNC(2), GF(2)
     
     
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
			 DNC=Norm2(cellCenter(IN,JN,:)-CellCenter(I,J,:)) 
             NF(:)=SF(IFace,:)/Norm2(SF(IFace,:))
			 
!for skew correction
             RNC(:)=(CellCenter(IN,JN,:)-CellCenter(I,J,:)) /DNC
			 GF(1)=RLinearInterp(DC,DN,GradP(I,J,1),GradP(IN,JN,1))
			 GF(2)=RLinearInterp(DC,DN,GradP(I,J,2),GradP(IN,JN,2))
			 dpdn=(P(IN,JN)-P(I,J))/DNC
			 
			 IF (DN.lt.1e-5) THEN
			    dpdn_c=DOT_PRODUCT(GradP(I,J,:),NF(:))
				dpdn=5./3*dpdn-2./3*dpdn_c
				GF(:)=GradP(I,J,:)
			 END IF
!skew correction
             dpdn=dpdn+DOT_PRODUCT(NF(:)-RNC(:),GF(:))
             
			 LapP(I,J)=LapP(I,J)+dpdn*Norm2(SF(IFace,:))
		  END DO
		  LapP(I,J)=LapP(I,J)/VOL
        END DO
     END DO
     
End Subroutine 
