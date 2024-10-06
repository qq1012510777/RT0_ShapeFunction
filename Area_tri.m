function f = Area_tri(p1, p2, p3)
    a = norm(p1 -p2);
    b = norm(p3 -p2);
    c = norm(p1 -p3);

    p = (a + b + c) / 2.0;

    f = (p * (p - a) * (p - b) * (p - c))^0.5;
end
