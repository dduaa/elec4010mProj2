%ELEC4010M Project2
%LI Junrong
%20328496

%prompt for inputs
prompt = 'What is the number of vertices? (>3) ';
n = input(prompt);
vertices = 1 : 1 : 2*n;

vertices_x = 1:1:n;
vertices_y = 1:1:n;

for i = 1 : n
    prompt = ['What is the X of vertex ', num2str(i),' ?: '];
    vertices(2*i-1) = input(prompt);
    vertices_x(i) = vertices(2*i-1);
    prompt = ['What is the Y of vertex ', num2str(i),' ?: '];
    vertices(2*i) = input(prompt);
    vertices_y(i) = vertices(2*i);
end
prompt = 'What is the friction coefficient ?: ';
%turn to angle
u = input(prompt);
u_angle = atan(u) * 180 / pi;

%number of line segments
s = 10;

%visualize the polygon
figure;
plot(vertices_x,vertices_y,'color','red');
hold on;
makeup_vertices_x = [vertices_x(1) vertices_x(n)];
makeup_vertices_y = [vertices_y(1) vertices_y(n)];
plot(makeup_vertices_x,makeup_vertices_y,'color','red');
title('Polygon visualization');
axis([min(vertices_x)-2 max(vertices_x)+2 min(vertices_y)-2 max(vertices_y)+2]);

% 2 contacts situation
% for i = 1 : n
%     for j = i+1 : n
%         %edge 1
%         x_1 = vertices_x(i);
%         y_1 = vertices_y(i);
%         
%         temp_index = i + 1;
%         if(temp_index > n)
%             temp_index = 1;
%         end
%         
%         x_2 = vertices_x(temp_index);
%         y_2 = vertices_y(temp_index);
%         
%         %edge 2
%         x_3 = vertices_x(j);
%         y_3 = vertices_y(j);
%         
%         temp_index = j + 1;
%         if(temp_index > n)
%             temp_index = 1;
%         end
%         x_4 = vertices_x(temp_index);
%         y_4 = vertices_y(temp_index);
%         
%         %compute the direction perpendicular to the line 
%         direction_1 = -1/((y_2-y_1)/(x_2-x_1));
%         angle_1 = atan(direction_1) * 180 / pi;
%         direction_2 = -1/((y_4-y_3)/(x_4-x_3));
%         angle_2 = atan(direction_2) * 180 / pi;
%         
%         plot_x = zeros(1,s);
%         plot_y = zeros(1,s);
%         plot_position_index = 1;
%         
%         for a = 1:1:s-1
%             %compute target point 1
%             t1_x = x_1 * (a/s) + x_2 * (1 - a/s);
%             t1_y = y_1 * (a/s) + y_2 * (1 - a/s);
%             
%             for b = 1:1:s-1 
%                 %compute target point 2
%                 t2_x = x_3 * (b/s) + x_4 * (1 - b/s);
%                 t2_y = y_3 * (b/s) + y_4 * (1 - b/s);
%                 connection_tan = (t2_y - t1_y)/(t2_x - t1_x);
%                 connection_angle = atan(connection_tan) * 180 / pi;
%                 compare_1 = abs(connection_angle-angle_1);
%                 compare_2 = abs(connection_angle-angle_2);
%                 if(and(compare_1 < u_angle,compare_2 < u_angle))
%                     plot_x(plot_position_index) = a/s;
%                     plot_y(plot_position_index) = b/s;
%                     plot_position_index = plot_position_index + 1;
%                 end
%             end
%         end
%         
%         figure;
%         hold on;
%         axis([1/s 1 1/s 1]);
%         title('2 contacts situation');
%         label = ['This is axis ', num2str(i)];
%         xlabel(label);
%         label = ['This is axis ', num2str(j)];
%         ylabel(label);
%         scatter(plot_x,plot_y);
%         
%     end
% end

% 3 contacts situation
for i = 1 : n-2
    for j = i+1 : n-1
        for k = j+1 : n
            %edge 1
            x_1 = vertices_x(i);
            y_1 = vertices_y(i);
        
            temp_index = i + 1;
            x_2 = vertices_x(temp_index);
            y_2 = vertices_y(temp_index);
            
            %edge 2
            x_3 = vertices_x(j);
            y_3 = vertices_y(j);
        
            temp_index = j + 1;
            x_4 = vertices_x(temp_index);
            y_4 = vertices_y(temp_index);
            
            %edge 3
            x_5 = vertices_x(k);
            y_5 = vertices_y(k);
        
            temp_index = k + 1;
            if(temp_index == n+1)
                temp_index = 1;
            end
            x_6 = vertices_x(temp_index);
            y_6 = vertices_y(temp_index);
            
            %compute the direction perpendicular to the line 
            direction_1 = -1/((y_2-y_1)/(x_2-x_1));
            angle_1 = atan(direction_1) * 180 / pi;
            direction_2 = -1/((y_4-y_3)/(x_4-x_3));
            angle_2 = atan(direction_2) * 180 / pi;
            direction_3 = -1/((y_6-y_5)/(x_6-x_5));
            angle_3 = atan(direction_3) * 180 / pi;
            
            %prepare for plotting
            plot_x_3 = zeros(1,s);
            plot_y_3 = zeros(1,s);
            plot_z_3 = zeros(1,s);
            plot_position_index = 1;
            
            for a = 1:1:s-1
                %compute target point 1
                t1_x = x_1 * (a/s) + x_2 * (1 - a/s);
                t1_y = y_1 * (a/s) + y_2 * (1 - a/s);
                for b = 1:1:s-1
                    %compute target point 2
                    t2_x = x_3 * (b/s) + x_4 * (1 - b/s);
                    t2_y = y_3 * (b/s) + y_4 * (1 - b/s);
                    for c = 1:1:s-1
                        %compute target point 3
                        t3_x = x_5 * (c/s) + x_6 * (1 - c/s);
                        t3_y = y_5 * (c/s) + y_6 * (1 - c/s);
                        
                        %check force closure
                        
                        %check Z
                        %we take t1 as the origin
                        test_direction_12 = (t2_y-t1_y)/(t2_x-t1_x);
                        test_angle_12 = atan(test_direction_12) * 180 / pi;
                        test_direction_13 = (t3_y-t1_y)/(t3_x-t1_x);
                        test_angle_13 = atan(test_direction_13) * 180 / pi;
                        
                        test_angle_1 = angle_2 + u_angle;
                        test_angle_2 = angle_2 - u_angle;
                        test_angle_3 = angle_3 + u_angle;
                        test_angle_4 = angle_3 - u_angle;
                        
                        test_angle_5 = angle_1 + u_angle;
                        test_angle_6 = angle_1 - u_angle;
                        
                        lock_positive = 0;
                        lock_negative = 0;
                        
                        if(test_angle_12 < test_angle_1)
                            lock_positive = 1;
                        elseif(test_angle_12 > test_angle_1)
                            lock_negative = 1;
                        end
                        if(test_angle_12 < test_angle_2)
                            lock_positive = 1;
                        elseif(test_angle_12 > test_angle_2)
                            lock_negative = 1;
                        end
                        if(test_angle_13 < test_angle_3)
                            lock_positive = 1;
                        elseif(test_angle_13 > test_angle_3)
                            lock_negative = 1;
                        end
                        if(test_angle_13 < test_angle_4)
                            lock_positive = 1;
                        elseif(test_angle_13 > test_angle_4)
                            lock_negative = 1;
                        end
                        if(lock_positive+lock_negative > 1)
                            mz_pass = 1;
                        else
                            mz_pass = 0;
                        end
                        
                        %check fx & fy
                        %called fxfx_pass
                        lock_phase_1 = 0;
                        lock_phase_2 = 0;
                        lock_phase_3 = 0;
                        lock_phase_4 = 0;
%                         disp(test_angle_1);
%                         disp(test_angle_2);
%                         disp(test_angle_3);
%                         disp(test_angle_4);
                        [temp_1,temp_2,temp_3,temp_4] = check_fxfy(test_angle_1);
                        lock_phase_1 = lock_phase_1 + temp_1;
                        lock_phase_2 = lock_phase_2 + temp_2;
                        lock_phase_3 = lock_phase_3 + temp_3;
                        lock_phase_4 = lock_phase_4 + temp_4;
                        [temp_1,temp_2,temp_3,temp_4] = check_fxfy(test_angle_2);
                        lock_phase_1 = lock_phase_1 + temp_1;
                        lock_phase_2 = lock_phase_2 + temp_2;
                        lock_phase_3 = lock_phase_3 + temp_3;
                        lock_phase_4 = lock_phase_4 + temp_4;
                        [temp_1,temp_2,temp_3,temp_4] = check_fxfy(test_angle_3);
                        lock_phase_1 = lock_phase_1 + temp_1;
                        lock_phase_2 = lock_phase_2 + temp_2;
                        lock_phase_3 = lock_phase_3 + temp_3;
                        lock_phase_4 = lock_phase_4 + temp_4;
                        [temp_1,temp_2,temp_3,temp_4] = check_fxfy(test_angle_4);
                        lock_phase_1 = lock_phase_1 + temp_1;
                        lock_phase_2 = lock_phase_2 + temp_2;
                        lock_phase_3 = lock_phase_3 + temp_3;
                        lock_phase_4 = lock_phase_4 + temp_4;
                        [temp_1,temp_2,temp_3,temp_4] = check_fxfy(test_angle_5);
                        lock_phase_1 = lock_phase_1 + temp_1;
                        lock_phase_2 = lock_phase_2 + temp_2;
                        lock_phase_3 = lock_phase_3 + temp_3;
                        lock_phase_4 = lock_phase_4 + temp_4;
                        [temp_1,temp_2,temp_3,temp_4] = check_fxfy(test_angle_6);
                        lock_phase_1 = lock_phase_1 + temp_1;
                        lock_phase_2 = lock_phase_2 + temp_2;
                        lock_phase_3 = lock_phase_3 + temp_3;
                        lock_phase_4 = lock_phase_4 + temp_4;
                        
                        fxfy_pass = 0;
                        if(lock_phase_1*lock_phase_2*lock_phase_3*lock_phase_4 > 0)
                            fxfy_pass = 1;
                        end
                        
                        if(and(mz_pass,fxfy_pass))
                            plot_x_3(plot_position_index) = a/s;
                            plot_y_3(plot_position_index) = b/s;
                            plot_z_3(plot_position_index) = c/s;
                            plot_position_index = plot_position_index + 1;
                        end
                        
                        
                    end
                end
            end
            
            figure;
            
            scatter3(plot_x_3,plot_y_3,plot_z_3);
            hold on;
            title('3 contacts situation');
            label = ['This is axis ', num2str(i)];
            xlabel(label);
            label = ['This is axis ', num2str(j)];
            ylabel(label);
            label = ['This is axis ', num2str(k)];
            zlabel(label);
            axis([1/s 1 1/s 1 1/s 1]);
            
        end
    end
end

