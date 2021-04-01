function result = bloodVesselConnection(seed,pix_di,lung,low_lung)
    Isizes = size(seed);
    bw_s = bwskel(seed);
    bw_endP = bwmorph3(bw_s,'endpoints');
    bw_cc = bwconncomp(seed);
    idx_end = find(bw_endP==1);
    endpoints_Idx = cell(1,1);
    pos = 0;
    bw_branch = bwmorph3(bw_s,'branchpoints');
    idx_br = find(bw_branch==1);
    
    for i = 1:size(idx_end,1)
        for k = 1 : bw_cc.NumObjects
            if ismember(idx_end(i),bw_cc.PixelIdxList{k})
                [x y z] = ind2sub([Isizes(1) Isizes(2) Isizes(3)],idx_end(i));
                endpoints_Idx{pos+1,1} = {[x y z],idx_end(i)};
                endpoints_Idx{pos+1,2} = k;
                pos = pos+1;
                break
            end
        end
    end
    for i = 1:size(idx_br,1)
        if ~ismember(idx_br(i),idx_end)
            for k = 1 : bw_cc.NumObjects
                if ismember(idx_br(i),bw_cc.PixelIdxList{k}) 
                    [x y z] = ind2sub([Isizes(1) Isizes(2) Isizes(3)],idx_br(i));
                    endpoints_Idx{pos+1,1} = {[x y z],idx_br(i)};
                    endpoints_Idx{pos+1,2} = k;
                    pos = pos+1;
                    break
                end
            end
        end
    end
    
    endpoints_Idx = sortrows(endpoints_Idx,2);
    neigb_pix = cell(size(endpoints_Idx,1),1);
    while pos > 0 
        idx = endpoints_Idx{pos};
        x = idx{1,1}(1);
        y = idx{1,1}(2);
        z = idx{1,1}(3);
        n_pos = 0;
        keep = {};
        for i = 1:size(endpoints_Idx,1)
            if pos ~= i
                r_idx = endpoints_Idx{i};
                xn = r_idx{1,1}(1);
                yn = r_idx{1,1}(2);
                zn = r_idx{1,1}(3);
                if zn < z -pix_di  || zn > z + pix_di
                    continue
                end
                cond = xn <= x+pix_di && xn >= x -pix_di && ...
                       yn <= y+pix_di && yn >= y -pix_di && ...
                       zn <= z+pix_di && zn >= z -pix_di ;
                if cond && r_idx{2}~=idx{2} 
                    keep{n_pos+1} = [xn yn zn];
                    n_pos = n_pos+1;
                end
            end
        end
        neigb_pix{pos,1} = idx{1,1};
        neigb_pix{pos,2} = keep;
        neigb_pix{pos,3} = length(keep);
        pos = pos -1;
    end
    neigb_pix = sortrows(neigb_pix,3);
    connect_point = zeros(Isizes)>0;
    marker = zeros(Isizes)>0;
    for i = 1 : size(neigb_pix,1)
        val = neigb_pix{i,1};
        nb = neigb_pix{i,2};
        if ~isempty(nb) 
            d = [];
            n_pos = 1;
            for j = 1 : size(nb,2)
                X = [val(1) val(2) val(3) ...
                    ;nb{j}(1) nb{j}(2) nb{j}(3)];
                d(n_pos,:) = [nb{j}(1) nb{j}(2) nb{j}(3) pdist(X,'euclidean')];
                n_pos = n_pos+1;
            end
            d = sortrows(d,4);
            for j = 1 : n_pos-1
                if marker(d(j,1),d(j,2),d(j,3)) == 0
                    x = val;
                    y = [d(j,1) d(j,2) d(j,3)];
                    if abs(val(3)-d(j,3))<=5
                        break
                    end
                    
                    [p_idx,miss_p] = connectP2P(x,y,Isizes,lung,low_lung);
                    
                    if miss_p < 0.8
                        break
                    end
                    se = strel('sphere',1);
                    p_idx = imdilate(p_idx,se)>0;
                    connect_point(p_idx)=1;
                    fprintf('\n%0.2f %% of Blood Vessel Connection',i/size(neigb_pix,1)*100);
                    marker(x(1),x(2),x(3)) = 1;
                    break
                end
            end
        end
    end

    result =  connect_point;

end


function [Im,val] = connectP2P(x,y,Isizes,lung,low_lung)
    Im= zeros(Isizes)>0;
    Im(x(1),x(2),x(3))=1;
    Im(y(1),y(2),y(3))=1;
    val = 0;
    h = 0;
    while ~isequal(x,y)
        a = y(1)-x(1);
        b = y(2)-x(2);
        c = y(3)-x(3);
        dx=rad2deg(atan2(a,b));
        dy=rad2deg(atan2(b,a));
        dz= rad2deg(atan2(c,a));
        x(1) = chDirection(dx,x(1));
        x(2) = chDirection(dy,x(2));
        x(3) = chDirection(dz,x(3));
        if lung(x(1),x(2),x(3))==1 && low_lung(x(1),x(2),x(3))==1 
            h = h+1;
        end
        Im(x(1),x(2),x(3))=1;
        val = val+1;
    end
    Im=Im>0;
    val = h/val;
end

function re = chDirection(da,val)
    if da < 0
        da = da + 360;
    end
    if da >= 30 && da <= 150
        val = val+1;
    elseif da >= 210 && da <= 330
        val = val -1;
    end
    re = val;
end