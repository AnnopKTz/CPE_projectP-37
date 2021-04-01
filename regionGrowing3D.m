function I_bw = regionGrowing3D(im,seed,mask)
    t_low = 250;
    I_bw = zeros(size(im{1},1),size(im{1},2),size(im,2));
    
    %convert variable
    im = groupCell(im);
    
    Isizes = size(I_bw);
    disp(Isizes);
    
    %neigb
    neigb3D = [-1 -1 -1; -1 0 -1 ; -1 1 -1 ; 0 -1 -1; 0 0 -1; 0 1 -1; 1 -1 -1; 1 0 -1;...
        1 1 -1;-1 -1 0; -1 0 0 ; -1 1 0 ; 0 -1 0; 0 1 0; 1 -1 0; 1 0 0; 1 1 0;-1 -1 1;...
        -1 0  1; -1 1 1 ; 0 -1 1; 0 0 1; 0 1 1; 1 -1 1; 1 0 1; 1 1 1];
    

    quere = zeros(1,3);
    pos = 0;
    I_bw = seed>0;
    seed_edge = edge3(seed,'approxcanny',1);
    
    for n = 1 : Isizes(3)
        for i = 1 : Isizes(1)
            for j = 1: Isizes(2)
                if seed_edge(i,j,n) == 1
                    quere(pos+1,:)=[i j n];
                    pos = pos+1;
                end
            end
        end
    end
    
    while pos >0
        x = quere(pos,1);
        y = quere(pos,2);
        z = quere(pos,3);
        crr_pos = pos;
        for i = 1 :size(neigb3D,1)
            nb = neigb3D(i,:);
            xn = x + nb(1);
            yn = y + nb(2);
            zn = z + nb(3);
            
            cond = zn<=Isizes(3)&& zn>=1 && xn<=Isizes(1) && yn<=Isizes(2) && xn>=1 && yn>= 1;
                
            if cond && im(xn,yn,zn)>=t_low && mask(xn,yn,zn)==1 && I_bw(xn,yn,zn) ~= 1 
                quere(pos+1,:)=[i j n];
                pos = pos+1;
                I_bw(i,j,n) = 1;
            end
        end
        quere(crr_pos,:) = [];
        pos = pos -1;
    end
    I_bw= I_bw>0;
    %check centroid for delete lung's issues
    status = regionprops3(I_bw,'Centroid','VoxelIdxList','PrincipalAxisLength','Image');
    q_remove = find(status.Centroid(:,1)<= Isizes(1)/2+25 & status.Centroid(:,1)>= Isizes(1)/2-25 );
    for i = 1 : size(q_remove)
        if q_remove(i)>size(status,1)
            break
        elseif status.Centroid(q_remove(i),2)<= Isizes(2)/2
            I_bw(status.VoxelIdxList{q_remove(i)}) = 0;
        end
    end
    
end

function g = groupCell(im)
    g = zeros((size(im{1},1)),(size(im{1},2)),size(im,2));
    for i = 1 : size(im,2)
        g(:,:,i) = im{i};
    end
end

