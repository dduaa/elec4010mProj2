%{
        w = [0 0 0];
            data = [1 1 0;
                1 -1 0;
                1 2 0;
                1 1 1;
                1 -1 1];
            T = [1;-1;2;2;0];
            learning_rate=0.1;
            for iteration = 1:1000
                t=0;
                for i=1:5
                    o = sum (data(i,:).*w);
                    diff = T(i) - o;
                    w = w + learning_rate*(diff).*data(i,:);
                    if abs(diff)<1e-5
                        t=t+1;
                    end
                end

                end
            end
    
function [a,b,c,d]=check_fxfy(check_angle)
a=0;
b=0;
c=0;
d=0;
t = check_angle;
if(t == -180)
    c = 1;
elseif(t<-90)
    c = 1;
    d = 1;
elseif(t==-90)
    d=1;
elseif(t<0)
    a=1;
    d=1;
elseif(t==0)
    a=1;
elseif(t<90)
    a=1;
    b=1;
elseif(t==90)
    b=1;
elseif(t<180)
    b=1;
    c=1;
elseif(t==180)
    c=1;
end

    

end

%}
function matrix
A=zeros(5,5);


for i=1:10000
    
    temp = A;
    for j=1:5
        for k=1:5
            
            if j==1 && k==2
                A(1,2) = 10 + temp(5,2);
                continue;
            end
            if j==1 && k==4
                A(1,4) = 5+temp(3,4);
                continue;
            end
            total = 0;
            total = total + 0.25*cal(temp,j+1,k,j,k);
            total = total + 0.25*cal(temp,j,k+1,j,k);
            total = total + 0.25*cal(temp,j-1,k,j,k);
            total = total + 0.25*cal(temp,j,k-1,j,k);
            A(j,k) = total;
            
        end
    end
    
end

A


    function value = cal( A , m,n,j,k)
       if m==0
            value = A(j,k)-1;
        elseif m==5
            value = A(j,k)-1;
        elseif n==0
            value = A(j,k)-1;
        elseif n==5
            value = A(j,k)-1;
        else
            value = A(m,n)-1;
        end
    end

end