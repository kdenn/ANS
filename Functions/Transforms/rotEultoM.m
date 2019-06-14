function R = rotEultoM(eul)
R = R3(eul(1))*R1(eul(2))*R3(eul(3));
end