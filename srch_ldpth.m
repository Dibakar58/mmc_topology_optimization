function [strct,loadPth] = srch_ldpth(nAct,allPhiAct,Phimax,epsilon,eleNodesID,loadEle,fixEle)
strct = 0;  fixSet = [];     loadPth= [];   Frnt= [];	% initialization
Hmax = Heaviside(Phimax,0,epsilon);	denmax = sum(Hmax(eleNodesID), 2)/8;
allH = Heaviside(allPhiAct,0,epsilon);	allden = repmat(denmax,1,nAct);
for i = 1:nAct
    tempH = allH(:,i);
    allden(:,i) = sum(tempH(eleNodesID), 2)/8;	% density matrix of all active components
end
if min(denmax(loadEle(:)))>eps && max(denmax(fixEle(:)))>eps
    cnnt=sparse(nAct,nAct);	% connection matrix of components
    for i = 1:nAct
        for j = i+1:nAct
            if max(min(allden(:,[i,j]),[],2))>0
                cnnt(i,j) = 1;	cnnt(j,i) = 1;
            end
        end
        if  max(allden(loadEle,i))>0
            loadPth = unique([loadPth i]);
            Frnt =  unique(setdiff([Frnt find(cnnt(i,:)>0)],loadPth));  % search front
        end
        if max(allden(fixEle,i))>0
            fixSet = unique([fixSet i]);
        end
        while ~isempty(Frnt)
            loadPth =  sort([loadPth Frnt]); Temp = [];
            for  i =  Frnt
                Temp = [Temp find(cnnt(i,:)>0)];
            end
            Frnt = unique(setdiff(Temp,loadPth));
        end
        if ~isempty(intersect(loadPth,fixSet))
            strct = 1;
        end
    end
end
end