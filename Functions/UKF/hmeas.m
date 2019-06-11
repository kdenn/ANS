function z = hmeas(xACAF,JD,ast,rSat,R,cam)
xACI = rotVert(xACAF',ast.alp,ast.del,ast.W_0,ast.W_d,JD,false);
z = (transPt32(xACI,rSat,R,cam))';
end