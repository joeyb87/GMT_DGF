wts_bis = zeros(size(wts));
for ll = 1:lmax
    startindex = 2*(ll^2 - 1) + 1;
    endindex = 2*((ll+1)^2 - 1);
    disp(['l = ' num2str(ll) ', start: ' num2str(startindex) ', end: ' num2str(endindex)]); 
    for mm = -ll:ll
        if mm < 0
            index2 = endindex - 2*(ll-abs(mm));
            index1 = index2 - 1;
            oldidx1 = startindex + 2*abs(abs(mm)-ll);
            oldidx2 = 1 + startindex + 2*abs(abs(mm)-ll);
            disp(['    m = ' num2str(mm) ', index = ' num2str(index1) ' and ' num2str(index2)])
            wts_bis(index1:index2) = wts(oldidx1:oldidx2);
        elseif mm > 0
            index1 = startindex + 2*abs(mm-ll);
            index2 = index1 + 1;
            oldidx1 = endindex - 2*(ll-abs(mm)) - 1;
            oldidx2 = endindex - 2*(ll-abs(mm));
            disp(['    m = ' num2str(mm) ', index = ' num2str(index1) ' and ' num2str(index2)])
            wts_bis(index1:index2) = wts(oldidx1:oldidx2);
        elseif mm == 0
            index1 = startindex + ll*2;
            index2 = index1 + 1;
            disp(['    m = ' num2str(mm) ', index = ' num2str(index1) ' and ' num2str(index2)])
            wts_bis(index1:index2) = wts(index1:index2);
        end
    end
end
figure; plot(abs(wts_bis))
        