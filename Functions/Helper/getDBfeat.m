function info = getDBfeat(ftDB,ftNum,tag)
%[x,y,z,JD,count,Ptri,Mx,My,sum_r,num,sum_alpha]

switch tag
    case 'posACAF'
        info = ftDB(ftNum,1:3)';
    case 'JD'
        info = ftDB(ftNum,4);
    case 'count'
        info = ftDB(ftNum,5);
    case '3Dcov'
        Pvec = ftDB(ftNum,6:11);
        info = [Pvec(1) Pvec(2)  Pvec(3);
                Pvec(2) Pvec(4)  Pvec(5);
                Pvec(3) Pvec(5) Pvec(6)];
    case '2Dcov'
        info = [ftDB(ftNum,12:13);
                ftDB(ftNum,13:14)];
    case 'r'
        info = ftDB(ftNum,15)/ftDB(ftNum,5);
    case 'num'
        info = ftDB(ftNum,16);
    case 'alpha'
        info = ftDB(ftNum,17);
    otherwise
        info = null;
end
end