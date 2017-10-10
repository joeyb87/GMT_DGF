function newordering = reverse_m_index(oldordering,lmax,whichcurrents,idxsize)

newordering = zeros(size(oldordering));

if idxsize == 1
    for ll = 1:lmax
        startindex = (ll^2 - 1) + 1;
        endindex = ((ll+1)^2 - 1);
        %                 disp(['l = ' num2str(ll) ', start: ' num2str(startindex) ', end: ' num2str(endindex)]);
        for mm = -ll:ll
            if mm < 0
                index2 = endindex - (ll-abs(mm));
                oldidx1 = startindex + abs(abs(mm)-ll);
                %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index2)])
                newordering(index2) = oldordering(oldidx1);
            elseif mm > 0
                index1 = startindex + abs(mm-ll);
                oldidx2 = endindex - (ll-abs(mm));
                %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index1)])
                newordering(index1) = oldordering(oldidx2);
            elseif mm == 0
                index1 = startindex + ll;
                %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index1)])
                newordering(index1) = oldordering(index1);
            end
        end
    end
else
    for ll = 1:lmax
        startindex = length(whichcurrents)*(ll^2 - 1) + 1;
        endindex = length(whichcurrents)*((ll+1)^2 - 1);
        %                 disp(['l = ' num2str(ll) ', start: ' num2str(startindex) ', end: ' num2str(endindex)]);
        for mm = -ll:ll
            if mm < 0
                index2 = endindex - 2*(ll-abs(mm));
                index1 = index2 - 1;
                oldidx1 = startindex + 2*abs(abs(mm)-ll);
                oldidx2 = 1 + startindex + 2*abs(abs(mm)-ll);
                %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index1) ' and ' num2str(index2)])
                newordering(index1:index2) = oldordering(oldidx1:oldidx2);
            elseif mm > 0
                index1 = startindex + 2*abs(mm-ll);
                index2 = index1 + 1;
                oldidx1 = endindex - 2*(ll-abs(mm)) - 1;
                oldidx2 = endindex - 2*(ll-abs(mm));
                %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index1) ' and ' num2str(index2)])
                newordering(index1:index2) = oldordering(oldidx1:oldidx2);
            elseif mm == 0
                index1 = startindex + ll*2;
                index2 = index1 + 1;
                %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index1) ' and ' num2str(index2)])
                newordering(index1:index2) = oldordering(index1:index2);
            end
        end
    end
    
end
