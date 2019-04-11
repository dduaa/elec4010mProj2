%{
c1=[1;1];
c2=[1;0];
c3 = [2;0];
f1 = [0;-1];
f2=[1;0];
f3=[-1;-1];
f4=[1;1];
f5=[-1;1];
f6=[1;1];
                            w1 = [0;f1(1);f1(2)];
                            w2 = [0;f2(1);f2(2)];
                            d2 = c2-c1;
                            m2 = cross(  [d2(1) d2(2) 0], [f1(1)  f1(2) 0] );
                            m22 = cross(  [d2(1) d2(2) 0] , [f2(1)  f2(2) 0] );
                            w3 = [ m2(3);f3(1);f3(2)];
                            w4 = [m22(3);f4(1);f4(2)];
                            d2 = c3 - c1;
                            m3 = cross(  [d2(1) d2(2) 0] , [f5(1)  f5(2) 0] );
                            m33 = cross(  [d2(1) d2(2) 0], [f6(1)  f6(2) 0] );
                            w5 = [m3(3);f5(1);f5(2)];
                            w6 = [m33(3);f6(1);f6(2)];
                            wrench = [w1 w2 w3 w4 w5 w6];
f = [1 1 1 1];
A = [[-1 0 0 0];
    [0 -1 0 0];
    [0 0 -1 0];
    [0 0 0 -1]];
b= [ -1 -1 -1 -1];
w=[w1 w2 w3 w4];
 [x,fval,exitflag,output] = linprog(f,A,b,w,zeros(3,1));
syms f  u1 u2 v1 v2 delta a1 a2 a3 dx dy dz;

eqn1 =  f*a1/a3 ==u1;
eqn2 = f*a2/a3 ==v1;
eqn3 = f*(a1+delta*dx)/(a3+delta*dz) == u2;
eqn4 = f*(a2+delta*dy)/(a3+delta*dz) ==v2;
eqn = [eqn1,eqn2, eqn3, eqn4];
[solf, sola1,sola2,sola3] = solve( eqn, [f a1 a2 a3]);

vertices = [0 1 2 3 3 3 3 2 1;
                  0 0 0 0 1 2 3 2 1];
    xq=[0 0 1];
    yq=[0 1 0];
[in] = inpolygon(xq,yq,vertices(1,:),vertices(2,:))
            inside_y = yq(in)
%}

vertices = [0 1 2 3 3 3 3 2 1;
                  0 0 0 0 1 2 3 2 1];
    xq=[3+0.05 3-0.05];
    yq=[0.5 0.5];
[in,on] = inpolygon(xq,yq,vertices(1,:),vertices(2,:))
 inside_y = yq((in - on)>0);
            inside_x = xq((in - on)>0)
           inside_y = yq((in-on)>0)
