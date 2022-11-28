t_span = 0:0.04: 1;

[t,y] = ode45(@(t,y)right_part(t,y), t_span, [0; 0; 2.28]);

plot(t,y(:,1));
grid on;
dlmwrite('Graph_3god.csv', [t,y],';')

t_array = [];
for A = 0.2:0.2:1
    
for i = 1: length(t)
    if (y(i,1) <= A)&&(y(i + 1,1) > A) 
        t_ = (t(i+1)-t(i))/(y(i + 1,1) - y(i,1)) * (A - y(i,1)) + t(i);
    break
    end
end

t_array = [t_array; t_];
end

dlmwrite('uzli.csv', t_array, ';')
dlmwrite('f_values', y(:,2),';')



















function dy = right_part(t, y)
dy = [y(2); y(3); -1.814*y(2)^2];
end


