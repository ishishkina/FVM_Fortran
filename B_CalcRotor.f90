Subroutine B_CalcRotor(NI,NJ,V,rotV, CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)

INTEGER I,J, NI,NJ,IN,JN
REAL CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),RN(2)
REAL IFaceCenter(NI,NJ-1,2),IFaceVector(NI,NJ-1,2),JFaceCenter(NI-1,NJ,2),JFaceVector(NI-1,NJ,2)
REAL rotV(0:NI,0:NJ), V(0:NI,0:NJ,2), NCELL(4,2),RF(4,2),SF(4,2),VOL,RC(2),DC,DN,VF(2)
     
     
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
             
             rotV(I,J)=rotV(I,J)+(VF(2)*SF(IFACE,1)-VF(1)*SF(IFACE,2))
         end do
       
		  rotV(I,J)=rotV(I,J)/VOL
           END DO
        END DO
      
End Subroutine 
