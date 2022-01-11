function [ returnVar ] = isover( inVar )
switch inVar
    case true
        returnVar = msgbox('You are not done yet!');
        return
    case false
        returnVar = msgbox('You are STILL not done yet!');
        return
    otherwise
        returnVar = msgbox('FINALLY!!');
        exit
end
end
