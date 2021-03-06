function [bordermask,bridgeflag,vIdxOrig,bridgeEnds]=splitdeflections_4_bwboundaries_test_global(orderedset,bordermask,nucr)
bridgeflag=0; %returned as 1 if deflections are bridged
nucr=round(nucr/4)*4; %make sure nucr is a multiple of 4
curvature_threshold = 2; % 2 originally
perilength=size(orderedset,1);
bridgeEnds = {};
%%% detect deflection vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vIdx=getdeflections_noclose(orderedset,nucr); %returns boundary indices
vIdxOrig = vIdx;
%%% count vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vnum=length(vIdx);
if vnum<2 | vnum>100
    return; %if less than two vertices are detected, exit function
end
%%% calculate perimeter distance between adjacent vertices %%%%%%%%%%%%%%%%
% periIdx=vIdx;
% periIdxadj1=[periIdx(2:end);perilength+periIdx(1)];
% pairperi1=periIdxadj1-periIdx;

vertX = orderedset(vIdx,1);
vertY = orderedset(vIdx,2);
vpos = [vertX , vertY];
distMat = squareform(pdist([vertX,vertY]));
% Avoid long bridges
distMat(distMat > nucr * 3) = Inf;
perimMat = zeros(size(distMat));
for i = 1:length(vIdx)
    periIdx1=vIdx(i);
    for j = 1:length(vIdx)
        periIdx2=vIdx(j);
        temp1 = periIdx1 - periIdx2;
        temp2 = -temp1;
        temp1 = (temp1 < 0) * perilength + temp1;
        temp2 = (temp2 < 0) * perilength + temp2;
        perimMat(i , j) = min([temp1 temp2]);
    end
end
curvMat = perimMat ./ distMat;
curvMat = triu(curvMat , 1);
% Keep track of bridged vertices
bridged = [];
bridged_adj_mat = zeros(length(vIdx));
%%% pair and bridge vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while vnum >= 2
    [bestcurv , curvelidx] = max(curvMat(:));
    [bridge_i , bridge_j] = ind2sub(length(vIdx) , curvelidx);

    % If best curvature is too low, stop segmenting this object.
    if bestcurv < curvature_threshold
        break;
    end

    bridgeflag = 1;

    %%% Bridge the vertices defining the best curvature %%%%%%%%%%%%%%%%%%%
    [bx,by]=bridge_NR(vpos(bridge_i , :) , vpos(bridge_j , :));
    % [bx,by] can equal NaN if the vertices were the same position (occurs
    % when two objects share a single coordinate. If this is the case,
    % skip mapping the bridge and proceed with perimeter update and vertex
    % removal.
    bridgeEnds = [bridgeEnds {[bx(1) by(1);bx(end) by(end)]}];
    if ~isnan(bx)
        for bci=1:length(bx)
            %bridgeflag=1;
            bordermask(by(bci),bx(bci))=1;
        end
    end

    %%% assign new perimeter distances & remove old vertices %%%%%%%%%%%%%%
%     vIdx([bridge_i , bridge_j]) = [];
%     vnum = length(vIdx);

    vnum = vnum - 2;
    bridged = [bridged , bridge_i , bridge_j];
    bridged_adj_mat(bridge_i , bridge_j) = 1;
    bridged_adj_mat(bridge_j , bridge_i) = 1;
    
    if vnum < 2
        break;
    end

%     vertX = orderedset(vIdx,1);
%     vertY = orderedset(vIdx,2);
%     vpos = [vertX , vertY];
%     distMat = squareform(pdist([vertX,vertY]));
%     % Avoid long bridges
%     distMat(distMat > nucr * 2) = Inf;
    
    perimMat = zeros(size(distMat));

    % Calculate perimeter distance, crossing bridge when possible
    for i = 1:length(vIdx)
        periIdx1=vIdx(i);
        for j = i+1:length(vIdx)
            if ismember(i , bridged) || ismember(j , bridged)
                distMat(i , j) = Inf;
                perimMat(i , j) = 0;
                continue;
            end

            temp1 = 0;
            temp2 = 0;
            
            % Go in + direction
            counter_pre = i;
            counter = i;
            just_crossed = 0;
            while counter_pre ~= j
                
                if ismember(counter_pre , bridged) && just_crossed == 0
                    counter = find(bridged_adj_mat(counter , :));
                    temp1 = temp1 + pdist(vpos([counter , counter_pre] , :));
                    just_crossed = 1;
                elseif counter_pre == length(vIdx)
                    counter = 1;
                    temp1 = temp1 + vIdx(counter) - vIdx(counter_pre) + perilength;
                    just_crossed = 0;
                else
                    counter = counter + 1;
                    temp1 = temp1 + vIdx(counter) - vIdx(counter_pre);
                    just_crossed = 0;
                end
                % If comes back to i 
                if counter == i
                    temp1 = 0;
                    break;
                end
                counter_pre = counter;                           
            end

            % Go in - direction
            counter_pre = i;
            counter = i;
            just_crossed = 0;
            while counter_pre ~= j
                
                if ismember(counter_pre , bridged) && just_crossed == 0
                    counter = find(bridged_adj_mat(counter , :));
                    temp2 = temp2 + pdist(vpos([counter , counter_pre] , :));
                    just_crossed = 1;
                elseif counter_pre == 1
                    counter = length(vIdx);
                    temp2 = temp2 + vIdx(counter_pre) - vIdx(counter) + perilength;
                    just_crossed = 0;
                else
                    counter = counter - 1;
                    temp2 = temp2 + vIdx(counter_pre) - vIdx(counter);
                    just_crossed = 0;
                end
                % If comes back to i 
                if counter == i
                    temp2 = 0;
                    break;
                end
                counter_pre = counter;
            end

            perimMat(i , j) = min([temp1 temp2]);
        end
    end
    curvMat = perimMat ./ distMat;
    curvMat = triu(curvMat , 1);

end

%{
while vnum>=2
    vpos=orderedset(vIdx,:);
    %%% Determine adjacent vertices that define the highest curvature %%%%%
    vposadj1=[vpos(2:end,:);vpos(1,:)];
    pair1=vposadj1-vpos;
    pairdist1=sqrt(sum(pair1.^2,2));
    curvature1=pairperi1./pairdist1;
    [bestcurve,curve1idx]=sort(curvature1);
    % If best curvature is too low, stop segmenting this object.
    if bestcurve(end)<2
        break
    end
    bestcurveidx=curve1idx(end);
    if bestcurveidx==vnum
        bestcurveidxadj=1;
    else
        bestcurveidxadj=bestcurveidx+1;
    end
    % If this point is reached, a split will be performed, so mark it.
    bridgeflag=1;
    %%% Bridge the vertices defining the best curvature %%%%%%%%%%%%%%%%%%%
    [bx,by]=bridge(vpos(bestcurveidx,:),vpos(bestcurveidxadj,:));
    % [bx,by] can equal NaN if the vertices were the same position (occurs
    % when two objects share a single coordinate. If this is the case,
    % skip mapping the bridge and proceed with perimeter update and vertex
    % removal.
    bridgeEnds = [bridgeEnds {[bx(1) by(1);bx(end) by(end)]}];
    if ~isnan(bx)
        for bci=1:length(bx)
            %bridgeflag=1;
            bordermask(by(bci),bx(bci))=1;
        end
    end

    %%% assign new perimeter distances & remove old vertices %%%%%%%%%%%%%%
    previdx=bestcurveidx-1;
    if previdx==0
        previdx=vnum;
    end
    % Given vertex 3-4 gave best curvature and is now bridged, define the
    % perimeter from vertex 2 to 5: p(2-5)=p(2-3)+bridge+p(4-5).
    pairperi1(previdx)=pairperi1(previdx)+length(bx)-1+pairperi1(bestcurveidxadj);
    % Remove the vertices and perimeters of the vertices defining the best
    % curve.
    vIdx([bestcurveidx,bestcurveidxadj])=[];
    pairperi1([bestcurveidx,bestcurveidxadj])=[];
    vnum=length(vIdx);
end
%}

%keyboard;
end
%%% debug: visualize deflections on boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
vpos=orderedset(vIdx,:);
for vc=1:size(vpos,1)
    bordermask(vpos(vc,2),vpos(vc,1))=1;
end

tempbridgeimage=zeros(size(scmask));
for bci=1:length(bcx)
   tempbridgeimage(bcy(bci),bcx(bci))=1;
end
tempimage=mat2gray(scmask);
tempimage(:,:,2)=tempbridgeimage;
tempimage(:,:,3)=0;
imshow(tempimage);
%}
