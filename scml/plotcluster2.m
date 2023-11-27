function [] = plotcluster2(X, lab)  
% This function plots the embedding.
%
%   Parameters are:
%
%   'X'     - N by D matrix. Each row in X represents an observation.
%   'lab'   - True annotations of data X.

if min(lab==fix(lab))==1 && all(lab >= 0)
    colors = [[31,119,179]/255;[251,130,20]/255;[43,159,46]/255;[210,33,33]/255;...
    [143,99,187]/255;[140,87,76]/255;[255,116,192]/255;[125,125,125]/255;...
    [184,187,29]/255;[30,191,208]/255;[218,165,32]/255;[65,105,225]/255;...
    [255,99,71]/255;[147,112,219]/255;[255,215,0]/255;[50,205,50]/255;...
    [174,199,232]/255;[255,187,120]/255;[152,223,138]/255;[255,152,150]/255;...
    [196,177,213]/255;[196,155,147]/255;[219,219,141]/255;[135,206,235]/255;...
    [255,165,0]/255;[144,238,144]/255;[1,0,0];[0,0,1];[0,1,0];[1,1,0];...
    [1,0,1];[0,1,1];[160,0,160]/255;[12,128,144]/255;[255,69,0]/255;...
    [140,86,75]/255;[160,82,45]/255;[0,139,139]/255;[175,238,238]/255;...
    [233,150,122]/255;[143,188,143]/255;[106,90,205]/255;[60,179,113]/255;...
    [220,20,60]/255;[65,105,225]/255;[147,112,219]/255;[20,206,209]/255];
    for i=1:length(lab)
        if(lab(i)==0)
            plot(X(i,1),X(i,2),'o','markerfacecolor',[180,180,180]/255,'markeredgecolor',[180,180,180]/255,'markersize',3);
            hold on;
        elseif(lab(i) <= 40)
            plot(X(i,1),X(i,2),'o','markerfacecolor',colors(lab(i),:),'markeredgecolor',colors(lab(i),:),'markersize',3);
            hold on;
        elseif(lab(i) > 40)
            plot(X(i,1),X(i,2),'o','markerfacecolor',colors(41+mod(lab(i)-41,7),:),'markeredgecolor',colors(41+mod(lab(i)-41,7),:),'markersize',3);
            hold on;
        end
    end
else
    disp('WARNING: clustering annotation must be a non-negative integer!');
end
end
