function [ names ] = output_names_multiIF_puncta_fixed( sigs,localbg,ring, puncta, thresh )
names={1,2,3 ;'x','y','nuclear area'};
ind=3;

for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        ind=ind+1;
        names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} 'mean']};
    end
end
for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        ind=ind+1;
        names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} 'median']};
    end
end
for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        if localbg{i}(j)
            ind=ind+1;
            names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} 'block bg']};
        end
    end
end
for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        if ring{i}(j)
            ind=ind+1;
            names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} 'cyto ring']};
        end
    end
end
for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        if puncta{i}(j)
            for t=1:length(thresh{i}{j})
                ind=ind+1;
                names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} 'puncta area_' num2str(thresh{i}{j}(t))]};
            end
        end
    end
end
for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        if puncta{i}(j)
            for t=1:length(thresh{i}{j})
                ind=ind+1;
                names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} 'puncta intensity_' num2str(thresh{i}{j}(t))]};
            end
        end
    end
end
for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        if puncta{i}(j)
            for t=1:length(thresh{i}{j})
                ind=ind+1;
                names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} 'puncta number_' num2str(thresh{i}{j}(t))]};
            end
        end
    end
end


