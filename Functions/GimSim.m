function [sat,G1,G2,G3,G4]=GimSim(ax,Theangle ,DCM_bg_0, OMEGAS , sat ,G1,G2,G3,G4)
% title(('VSCMGs Control byAliAlmorshdy'))
delete(sat);delete(G1);delete(G2);delete(G3);delete(G4);
DCM_bg11      = DCM_bg_0(1:3,1:3);
DCM_bg22      = DCM_bg_0(4:6,1:3);
DCM_bg33      = DCM_bg_0(7:9,1:3);
DCM_bg44      = DCM_bg_0(10:12,1:3);
gs1=DCM_bg11(:,1);gt1=DCM_bg11(:,2);gg1=DCM_bg11(:,3);
DCM_bg1=[gt1,gg1,gs1];

gs2=DCM_bg22(:,1);gt2=DCM_bg22(:,2);gg2=DCM_bg22(:,3);
DCM_bg2=[gt2,gg2,gs2];

gs3=DCM_bg33(:,1);gt3=DCM_bg33(:,2);gg3=DCM_bg33(:,3);
DCM_bg3=[gt3,gg3,gs3];

gs4=DCM_bg44(:,1);gt4=DCM_bg44(:,2);gg4=DCM_bg44(:,3);
DCM_bg4=[gt4,gg4,gs4];

% Theangle=[45 45 45];
BN=ConvertAttitude(Theangle'*pi/180,'321','DCM');
q= quaternion(Theangle,"eulerd","ZYX","frame");
sat=poseplot(ax,q,[0 0 0],MeshFileName='AliSat.STL',ScaleFactor=1,PatchFaceColor='y',PatchFaceAlpha=0.2);
hold on
MatR1=rotz(OMEGAS(1))'*DCM_bg1'*BN;
MatR2=rotz(OMEGAS(1))'*DCM_bg2'*BN;
MatR3=rotz(OMEGAS(1))'*DCM_bg3'*BN;
MatR4=rotz(OMEGAS(1))'*DCM_bg4'*BN;
angles1=ConvertAttitude(MatR1,'DCM','321')'*180/pi;
angles2=ConvertAttitude(MatR2,'DCM','321')'*180/pi;
angles3=ConvertAttitude(MatR3,'DCM','321')'*180/pi;
angles4=ConvertAttitude(MatR4,'DCM','321')'*180/pi;
q1= quaternion([angles1(1) angles1(2) angles1(3)],"eulerd","ZYX","frame");
q2= quaternion([angles2(1) angles2(2) angles2(3)],"eulerd","ZYX","frame");
q3= quaternion([angles3(1) angles3(2) angles3(3)],"eulerd","ZYX","frame");
q4= quaternion([angles4(1) angles4(2) angles4(3)],"eulerd","ZYX","frame");
G1=poseplot(ax,q1,BN'*[.17 .17 -.17]',MeshFileName='ARW.STL',ScaleFactor=2.5,PatchFaceColor='r',PatchFaceAlpha=1);
 G2=poseplot(ax,q2,BN'*[.17 -.17 -.17]',MeshFileName='ARW.STL',ScaleFactor=2.5,PatchFaceColor='r',PatchFaceAlpha=1);
G3=poseplot(ax,q3,BN'*[-.17 .17 -.17]',MeshFileName='ARW.STL',ScaleFactor=2.5,PatchFaceColor='r',PatchFaceAlpha=1);
 G4=poseplot(ax,q4,BN'*[-.17 -.17 -.17]',MeshFileName='ARW.STL',ScaleFactor=2.5,PatchFaceColor='r',PatchFaceAlpha=1);
end