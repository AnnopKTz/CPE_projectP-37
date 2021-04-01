clear;

input_dir = "..\FUMPE\CT_scans\PAT018";
output_dir = "..\PAT018";
name_file = "OUT_18";


%_________________________MAIN__________________________
k = dir(input_dir);
keep = readData(k);

%450
lung = segLung(keep,450);
writeData(lung,"..\preLung_18",name_file,"png","arr")

%400
out = cutTh(keep,lung,400);
out = groupCell(out);

result = regionGrowing3D(keep,out,lung);
low_lung = cutLung(keep,220);
final = bloodVesselConnection(result,10,lung,low_lung);
x = result;
result(final)=1;
final_result = bwareaopen(result,50);
writeData(x,"..\pre_connect_18",name_file,"png","arr")
writeData(final_result,output_dir,name_file,"png","arr")



% % _________________________measure
% gt_dir = "..\..\GT18\dicom2\ScalarVolume_24\";
% k = dir(gt_dir);
% k = readData(k);
% out = groupCell(k)>0;
% figure(1);volshow(final_result);
% figure(2);volshow(out);
% ans = dice(final_result,out)




%__________________________FUNCTION______________________

function mask = segLung(k,th)
    mask=cell(1,size(k,2));
    for i = 1 : size(mask,2)
        mask{i} = k{i}<=th;
    end
    for i = 1 : size(mask,2)
        test = imfill(~mask{i},'holes');
        mask{i} = mask{i}- ~(test);
        if i <= size(k,2)*0.25
            %close bay
            se = strel('line',6,90);
            mask{i} = imclose(mask{i},se);
            se = strel('line',6,0);
            mask{i} = imclose(mask{i},se);
            %remove small area
            bw_out = bwareaopen(~mask{i},400);
            mask{i} = imfill(mask{i},'holes');
            mask{i} = bitand(mask{i},~bw_out);
            
        else
            mask{i} = bwareaopen(mask{i},150);
            cc = bwconncomp( mask{i});
            obj = cc.NumObjects;
            while obj>=2
                numPixels = cellfun(@numel,cc.PixelIdxList);
                [biggest,idx1] = max(numPixels);
                [smallest,idx2] = min(numPixels);
                %remove noise and Trachea
                if obj == 2 
                    if biggest*0.3 >= smallest
                        mask{i}(cc.PixelIdxList{idx2})= 0;
                    end
                else
                    mask{i}(cc.PixelIdxList{idx2})= 0;
                    cc = bwconncomp( mask{i});
                end
                obj = obj-1;
            end
            se = strel('line',15,0);
            mask{i} = imclose(mask{i},se);
            se = strel('line',25,110);
            mask{i} = imclose(mask{i},se);
            se = strel('line',15,70);
            mask{i} = imclose(mask{i},se);
            mask{i} = imfill(mask{i},'holes');
        end
    end
    
    mask = groupCell(mask);
    mask = mask>0;
    CC = bwconncomp(mask);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx1] = max(numPixels);
    numPixels(idx1) = 0;
    [big,idx2] = max(numPixels);
    for i = 1 : CC.NumObjects
        if i ~= idx1 && i ~= idx2
            mask(CC.PixelIdxList{i})=false;
        elseif i == idx2 && biggest*0.3 >= big
            mask(CC.PixelIdxList{i})=false;
        end
    end
end

function keep = readData(k)
    keep=cell(1,size(k,1)-2);
    for i = 1 : size(k,1)
        if i >= 3
            t = strcat(k(i).folder,'\',k(i).name);
            keep{i-2} =dicomread(dicominfo(t));
        end
    end
end

function writeData(k,dir,name,lname,typ)
    if(exist(dir)==0)
        mkdir(dir)
    end
    if typ == "c"
        for i = 1 : size(k,2)
            if i < 10
                imwrite(k{i},sprintf('%s/%s_000%d.%s',dir,name,i,lname));
            elseif i <100
                imwrite(k{i},sprintf('%s/%s_00%d.%s',dir,name,i,lname));
            else
                imwrite(k{i},sprintf('%s/%s_0%d.%s',dir,name,i,lname));
            end
        end
    else
        for i = 1 : size(k,3)
            if i < 10
                imwrite(k(:,:,i),sprintf('%s/%s_000%d.%s',dir,name,i,lname));
            elseif i <100
                imwrite(k(:,:,i),sprintf('%s/%s_00%d.%s',dir,name,i,lname));
            else
                imwrite(k(:,:,i),sprintf('%s/%s_0%d.%s',dir,name,i,lname));
            end
        end
    end
end

function out = cutTh(im,mask,th)
    out = cell(1,size(im,2));
    se = strel('disk',3);
    for i = 1 : size(im,2)
        rm_edge = imerode(mask(:,:,i),se);
        out{i} = im{i} >= th;
        out{i} = bitand(out{i},rm_edge);
    end
end

function g = groupCell(im)
    g = zeros((size(im{1},1)),(size(im{1},2)),size(im,2));
    for i = 1 : size(im,2)
        g(:,:,i) = im{i};
    end
end

function mask = cutLung(k,th)
    mask=cell(1,size(k,2));
    for i = 1 : size(mask,2)
        mask{i} = k{i}>=th;
    end
    mask = groupCell(mask);
    mask = mask>0;
end
