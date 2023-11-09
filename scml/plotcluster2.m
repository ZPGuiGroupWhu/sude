function [  ] = plotcluster2(X,lab)  
% This function plots the embedding.
for i=1:length(lab)
    if(lab(i)==0||lab(i)==-1)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[180,180,180]/255,'markeredgecolor',[180,180,180]/255,'markersize',3);
        hold on;
    elseif(lab(i)==1)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[31,119,179]/255,'markeredgecolor',[31,119,179]/255,'markersize',3);
        hold on;
    elseif(lab(i)==2)
          plot(X(i,1),X(i,2),'o','markerfacecolor',[251,130,20]/255,'markeredgecolor',[251,130,20]/255,'markersize',3);
          hold on;
    elseif(lab(i)==3)
          plot(X(i,1),X(i,2),'o','markerfacecolor',[43,159,46]/255,'markeredgecolor',[43,159,46]/255,'markersize',3);
          hold on;
    elseif(lab(i)==4)
          plot(X(i,1),X(i,2),'o','markerfacecolor',[210,33,33]/255,'markeredgecolor',[210,33,33]/255,'markersize',3);
          hold on;
    elseif(lab(i)==5)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[143,99,187]/255,'markeredgecolor',[143,99,187]/255,'markersize',3);
        hold on;
    elseif(lab(i)==6)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[140,87,76]/255,'markeredgecolor',[140,87,76]/255,'markersize',3);
        hold on;
    elseif(lab(i)==7)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[255,116,192]/255,'markeredgecolor',[255,116,192]/255,'markersize',3);
        hold on;
    elseif(lab(i)==8)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[125,125,125]/255,'markeredgecolor',[125,125,125]/255,'markersize',3);
        hold on;
    elseif(lab(i)==9)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[184,187,29]/255,'markeredgecolor',[184,187,29]/255,'markersize',3);
        hold on;
     elseif(lab(i)==10)
         plot(X(i,1),X(i,2),'o','markerfacecolor',[30,191,208]/255,'markeredgecolor',[30,191,208]/255,'markersize',3);
         hold on;
    elseif(lab(i)==11)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[218,165,32]/255,'markeredgecolor',[218,165,32]/255,'markersize',3);
        hold on;
    elseif(lab(i)==12)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[65,105,225]/255,'markeredgecolor',[65,105,225]/255,'markersize',3);
        hold on;  
     elseif(lab(i)==13)
         plot(X(i,1),X(i,2),'o','markerfacecolor',[255,99,71]/255,'markeredgecolor',[255,99,71]/255,'markersize',3);
         hold on; 
    elseif(lab(i)==14)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[147,112,219]/255,'markeredgecolor',[147,112,219]/255,'markersize',3);
        hold on;
    elseif(lab(i)==15)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[255,215,0]/255,'markeredgecolor',[255,215,0]/255,'markersize',3);
        hold on;
    elseif(lab(i)==16)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[50,205,50]/255,'markeredgecolor',[50,205,50]/255,'markersize',3);
        hold on;  
    elseif(lab(i)==17)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[174,199,232]/255,'markeredgecolor',[174,199,232]/255,'markersize',3);
        hold on; 
    elseif(lab(i)==18)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[255,187,120]/255,'markeredgecolor',[255,187,120]/255,'markersize',3);
        hold on;
    elseif(lab(i)==19)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[152,223,138]/255,'markeredgecolor',[152,223,138]/255,'markersize',3);
        hold on;  
    elseif(lab(i)==20)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[255,152,150]/255,'markeredgecolor',[255,152,150]/255,'markersize',3);
        hold on; 
    elseif(lab(i)==21)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[196,177,213]/255,'markeredgecolor',[196,177,213]/255,'markersize',3);
        hold on;  
    elseif(lab(i)==22)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[196,155,147]/255,'markeredgecolor',[196,155,147]/255,'markersize',3);
        hold on; 
    elseif(lab(i)==23)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[219,219,141]/255,'markeredgecolor',[219,219,141]/255,'markersize',3);
        hold on;
    elseif(lab(i)==24)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[135,206,235]/255,'markeredgecolor',[135,206,235]/255,'markersize',3);
        hold on; 
    elseif(lab(i)==25)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[255,165,0]/255,'markeredgecolor',[255,165,0]/255,'markersize',3);
        hold on;  
    elseif(lab(i)==26)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[144,238,144]/255,'markeredgecolor',[144,238,144]/255,'markersize',3);
        hold on; 
    elseif(lab(i)==27)
        plot(X(i,1),X(i,2),'o','markerfacecolor','r','markeredgecolor','r','markersize',3);
        hold on; 
    elseif(lab(i)==28)
        plot(X(i,1),X(i,2),'o','markerfacecolor','b','markeredgecolor','b','markersize',3);
        hold on; 
    elseif(lab(i)==29)
        plot(X(i,1),X(i,2),'o','markerfacecolor','g','markeredgecolor','g','markersize',3);
        hold on; 
    elseif(lab(i)==30)
        plot(X(i,1),X(i,2),'o','markerfacecolor','y','markeredgecolor','y','markersize',3);
        hold on; 
    elseif(lab(i)==31)
        plot(X(i,1),X(i,2),'o','markerfacecolor','m','markeredgecolor','m','markersize',3);
        hold on; 
    elseif(lab(i)==32)
        plot(X(i,1),X(i,2),'o','markerfacecolor','c','markeredgecolor','c','markersize',3);
        hold on; 
    elseif(lab(i)==33)
          plot(X(i,1),X(i,2),'o','markerfacecolor',[160,0,160]/255,'markeredgecolor',[160,0,160]/255,'markersize',3);
          hold on;
    elseif(lab(i)==34)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[12,128,144]/255,'markeredgecolor',[12,128,144]/255,'markersize',3);
        hold on;
    elseif(lab(i)==35)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[255,69,0]/255,'markeredgecolor',[255,69,0]/255,'markersize',3);
        hold on;
    elseif(lab(i)==36)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[140,86,75]/255,'markeredgecolor',[140,86,75]/255,'markersize',3);
        hold on;
    elseif(lab(i)==37)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[160,82,45]/255,'markeredgecolor',[160,82,45]/255,'markersize',3);
        hold on;
    elseif(lab(i)==38)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[0,139,139]/255,'markeredgecolor',[0,139,139]/255,'markersize',3);
        hold on;
     elseif(lab(i)==39)
         plot(X(i,1),X(i,2),'o','markerfacecolor',[175,238,238]/255,'markeredgecolor',[175,238,238]/255,'markersize',3);
         hold on;
    elseif(lab(i)==40)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[233,150,122]/255,'markeredgecolor',[233,150,122]/255,'markersize',3);
        hold on;
    elseif(lab(i)==41)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[143,188,143]/255,'markeredgecolor',[143,188,143]/255,'markersize',3);
        hold on;  
     elseif(lab(i)==42)
         plot(X(i,1),X(i,2),'o','markerfacecolor',[106,90,205]/255,'markeredgecolor',[106,90,205]/255,'markersize',3);
         hold on; 
    elseif(lab(i)>43&&mod(lab(i),6)==0)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[60,179,113]/255,'markeredgecolor',[60,179,113]/255,'markersize',3);
        hold on;
    elseif(lab(i)>43&&mod(lab(i),6)==1)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[220,20,60]/255,'markeredgecolor',[220,20,60]/255,'markersize',3);
        hold on;
    elseif(lab(i)>43&&mod(lab(i),6)==2)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[65,105,225]/255,'markeredgecolor',[65,105,225]/255,'markersize',3);
        hold on;
    elseif(lab(i)>43&&mod(lab(i),6)==3)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[147,112,219]/255,'markeredgecolor',[147,112,219]/255,'markersize',3);
        hold on;
    elseif(lab(i)>43&&mod(lab(i),6)==4)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[20,206,209]/255,'markeredgecolor',[20,206,209]/255,'markersize',3);
        hold on;
    elseif(lab(i)>43&&mod(lab(i),6)==5)
        plot(X(i,1),X(i,2),'o','markerfacecolor',[255,105,180]/255,'markeredgecolor',[255,105,180]/255,'markersize',3);
        hold on;
    end

end
end

