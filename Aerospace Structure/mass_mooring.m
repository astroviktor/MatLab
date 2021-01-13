%mass calculation mooring
%load the .mat file first!
function mt=mass_mooring(PD)
mt=0;
for i=1:length(NE)
    PD.ElmConnect(i,:)=PD.NodePos(i); 
    PD.ElmConnect(:,i)=PD.NodePos(1+1);
  if PD.NodePos(i,3)==PD.NodePos(i+1,3)
      if NodePos(i,2)==NodePos(i+1,2)
          %case in which the bar is x-axis orientated
          L(i)=abs(NodePos(i,1)-NodePos(1+1,1));
      elseif NodePos(i,1)==NodePos(i+1,1)
          %case in which the bar is y-axis orientated
          L(i)=abs(NodePos(i,2)-NodePos(i+1,2));
      end
  else
      %in this case the nodes are in a different z position
      %so we use the pitagora theorem to compute the length of the bar
      k(i)=NodePos(i,1)-NodePos(i+1,1);
      j(i)=NodePos(i,2)-NodePos(i+1,2);
      L(i)=sqrt((k(i)^2)+(j(i)^2));
  end
  if ElmMats(i)==1
      m(i)=MatSets(6).rho*MatSets(6).A*L(i);
  else
      m(i)=MatSets(1).rho*MatSets(1).A*L(i);
  end
      
mt=mt+m(i);
end


