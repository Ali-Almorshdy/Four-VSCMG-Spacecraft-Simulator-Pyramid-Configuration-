xx=axes();
eull=ConvertAttitude(log_sigmaBN,'MRP','321')';
title("VSCMG Conrol")
[sat,G1,G2,G3,G4]=GimSim(xx,eull(end,1:3)*180/pi,log_frames{1},wheelangles(end,1:4)*180/pi,"","","","","");
for i=1:100:length(eull)
[sat,G1,G2,G3,G4]=GimSim(xx,eull(i,1:3)*180/pi,log_frames{i},wheelangles(i,1:4)*180/pi,sat,G1,G2,G3,G4);
axis off
drawnow
end