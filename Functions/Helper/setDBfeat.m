function ftDB = setDBfeat(ftDB,ftNum,tag,info)
%[x,y,z,JD,count,Ptri,Mx,My,sum_r,num,sum_alpha]

switch tag
    case 'posACAF'
        ftDB(ftNum,1:3) = info';
    case 'JD'
        ftDB(ftNum,4) = info;
    case 'count'
        ftDB(ftNum,5) = info;
    case '3Dcov'
        ftDB(ftNum,6:11) = info([1:3,5,6,9]);                
    case '2Dcov'
        ftDB(ftNum,12:14) = info([1,2,4]);
    case 'r'
        ftDB(ftNum,15) = ftDB(ftNum,15) + info;
    case 'num'
        ftDB(ftNum,16) = info;
    case 'alpha'
        ftDB(ftNum,17) = info;
    otherwise
end
end