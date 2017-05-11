function v=rootpath(cellid)

v = [];
while cellid > 0
    v = [v cellid];    
    cellid = floor(cellid/2);
end

v = sort(v);