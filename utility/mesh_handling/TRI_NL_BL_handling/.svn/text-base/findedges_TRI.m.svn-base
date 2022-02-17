

idx = find(inpolygon(x,y,roi(:,1),roi(:,2)));

xx = x(idx);
yy = y(idx);
zz = 0*x(idx)+0.00001*rand(size(x(idx)));

TRIx = delaunay(xx,yy);

disp('looking for edges')
idxedge = [];
for ll=1:length(xx)
    
    idxtemp = find(ismember(TRIx(:,1),ll)+ismember(TRIx(:,2),ll)+ismember(TRIx(:,3),ll));
    %length(idxtemp)
    thet=0;
    for dd=1:length(idxtemp)
        pair = setdiff(TRIx(idxtemp(dd),:),ll);
        veca = [xx(pair(1))-xx(ll); yy(pair(1))-yy(ll); zz(pair(1))-zz(ll)];
        veca = veca./sqrt(dot(veca,veca));
        vecb = [xx(pair(2))-xx(ll); yy(pair(2))-yy(ll); zz(pair(2))-zz(ll)];
        vecb = vecb./sqrt(dot(vecb,vecb));
        thet = thet + acos(abs(dot(veca,vecb)));
    end
    %thet*180/pi
    if thet<270/180*pi;
        idxedge = cat(1,idxedge,ll);
        %thet;
    end
end


idxb = idx(idxedge); 
%idxbt = []; 
idxin = setdiff(idx,idxb);




clear TRISin tempdiff idxtemp idxtemp2 idxedge  ztop zbottom condicio1 condicion2


%{

%}


clear TRISin tempdiff %idxtemp idxtemp2  ztop zbottom