load cyto.mat

N=size(drop,1);



for cont=1:N
    if ~isreal(drop(cont,:))
        disp(['the sample n. ' num2str(cont) ' contains complex values']);
    end
end