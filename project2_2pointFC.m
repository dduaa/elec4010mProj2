function project2_2pointFC
clear all;
close all;
clc;
%% MATLAB version: R2016a
%{
            Author: DU DONGHONG
            20328343
            dduaa@connect.ust.hk
            2018 SPRING ELEC4010M  PROJECT2
            extra interactive choosing point
%}

[vertices,segment, v_num, my_unit ,u] = basic;

disp('Please just drag the legend away if it blocks something');
disp('We do not consider the vertice point for its unstability');
disp('Red line for shape');
disp('Blue o for the contact you choose');
disp('Blue x for corresponding contact that could form a 2 contact points FC.')
disp('If there is no x, it means no solution for this contact');
contact = zeros(v_num, 2, segment-1);
x=[];   y = [];


for  i = 1: v_num
    xs = vertices(1,i);         xe = vertices(1,1+mod(i+1,v_num));
    ys = vertices(2,i);         ye = vertices(2,1+mod(i+1,v_num));
    x=[x xs];   y=[y ys];
    x_interval = (xe-xs) / segment;
    y_interval = (ye-ys) / segment;
    for j =1:segment -1
        contact( i , 1 , j) = xs + x_interval * j;
        contact( i , 2 , j) = ys + y_interval * j;
        x=[x  contact( i , 1 , j)];
        y=[y  contact( i , 2 , j)];
    end
    x=[x xe];
    y=[y ye];
end

%% basic of GUI
FigH = figure('position',[360 500 400 400],'name','find the possible solution for a contact for FC');
axes( 'units','pixels', 'position',[100 50 200 200]);
LineH = plot(x,y,'r');
title(['unit is ' my_unit]);
axis equal;
xlabel('X') % x-axis label
ylabel('Y') % y-axis label
legend('shape');
TextH_e = uicontrol('style','text','position',[170 310 60 15]);
TextH_e.String =sprintf('edge# 1');

TextH_c = uicontrol('style','text','position',[170 350 60 15]);
TextH_c.String =sprintf('contact# 1');

SliderC = uicontrol('style','slider','position',[100 330 200 20],'min', 1, 'max', segment-1,'value',1,'SliderStep',[1/(segment-1) 1], 'Callback', @callbackfnC);

SliderE = uicontrol('style','slider','position',[100 290 200 20], 'min', 1, 'max', v_num,'value',1,'SliderStep',[1/(v_num-1) 1], 'Callback', @callbackfnE);

if(~exist('quality1Button','var'))
    quality1Button = uicontrol('Style','radiobutton','String','quality measure extra','pos',[0 0 200 25],'parent',FigH,'Callback', @callbackfnQ1);
end
set(quality1Button,'Value',0);

Contact = 1;   Edge = 1;

movegui(FigH, 'center')
replot;
%% real-time response

    function callbackfnQ1(source, eventdata)
        if(get(quality1Button,'Value')==1)
            replot;
        else
            close(figure(3));
        end
    end
   
    function callbackfnE(source, eventdata)
        Edge =source.Value;
        Edge = round(Edge);
        set(source,'Value',  Edge);
        TextH_e.String =sprintf('edge# %d',Edge);
        replot;
    end

    function callbackfnC(source, eventdata)
        Contact =source.Value;
        Contact = round(Contact);
        set(source,'Value',  Contact);
        TextH_c.String =sprintf('contact# %d',Contact);
        replot;
    end

    function replot
        if(get(quality1Button,'Value')==1)
                fig3= figure(3);
                cla(fig3);
                set(fig3,'name','FC quality measure with 2 contacts, angle method');
        end

         LineH = figure(1);
         hold off;
        LineH = plot(x,y,'r');
        hold on;
        LineH = plot( contact(Edge ,1 ,Contact),contact(Edge ,2 ,Contact),'ob');
        title(['unit is ' my_unit]);
        
        axis equal;
        xlabel('X') % x-axis label
        ylabel('Y') % y-axis label
        legend('shape','chosen contact');
        %% find corresponding FC contact point
        max_q = 0;
        for m = 1:v_num
            if m~= Edge
                for p = 1: segment -1
                    x_point_one = contact( m, 1, p );                  y_point_one = contact( m, 2, p );
                    x_point_two = contact(Edge ,1 ,Contact);     y_point_two = contact(Edge ,2 ,Contact);
                    angle_i = angle (x_point_one, y_point_one,x_point_two,y_point_two, vertices(1,m),vertices(2,m) );
                    angle_j = angle (x_point_two, y_point_two,x_point_one,y_point_one, vertices(1,Edge),vertices(2,Edge) );
                    if  angle_i + 1e-5 <  atan(u) && angle_j + 1e-5 <  atan(u)
                        LineH = figure(1);
                        hold on;
                        LineH = plot( contact(m ,1 ,p),contact(m ,2 ,p),'xb');
                        legend('shape','chosen contact','solution');
                        if(get(quality1Button,'Value')==1)
                            %% angle percentage
                            distance1 = sqrt(         (x_point_two-vertices(1,Edge))^2+(y_point_two-vertices(1,Edge))^2  );
                            distance2 = sqrt(     (x_point_one-vertices(1,m))^2+(y_point_one-vertices(1,m))^2  );
                            fig3= figure(3);
                            hold on;
                            z=exp(-angle_j)+exp(-angle_i);
                            if z > max_q
                                max_q = z;
                            end
                            stem3(distance1,distance2,z);
                            view(3);
                        end
                    end
                end
            end
        end
            if(get(quality1Button,'Value')==1)
                fig3 = figure(3);
                view(3);
                xlabel('edge chosen first') % x-axis label
                ylabel('edge chosen second') % y-axis label
                fprintf('The maximum quality is %d .\n',max_q);
            end
        
    end

    function  Angle  = angle( xs , ys , xe , ye, x_edge , y_edge )
        Angle1 = atan2((ye-ys),(xe-xs));
        Angle2 = atan2((y_edge-ys),(x_edge-xs));
        Angle = abs(Angle1-Angle2);
        if Angle >= pi*1.5 + 1e-5
            Angle = Angle - pi/2-pi;
        elseif Angle >=pi + 1e-5
            Angle = pi*1.5-Angle;
        elseif Angle >=pi/2 + 1e-5
            Angle = Angle - pi/2;
        else
            Angle = pi/2 - Angle;
        end
    end
    function [vertices,segment, v_num, my_unit ,u] = basic
        %% basic requests for vertieces
        v_num  =   str2double( input('How many vertices? Please enter a valid integer','s'));
        while isnan(v_num) || fix(v_num)~= v_num || v_num<=2
            v_num  =   str2double( input('How many vertices? Please enter an integer','s'));
        end
        vertices = zeros (2,v_num);
        for i=1:v_num
            disp(['For #',num2str(i),'vertex']);
            vertices(1,i) = str2double( input(' x is','s')    );
            while isnan( vertices(1,i) )
                vertices(1,i) = str2double( input('Not a number, try again. x is','s')    );
            end
            
            vertices(2,i) = str2double( input('y is','s')    );
            while isnan(vertices(2,i) )
                
                vertices(2,i) = str2double( input('Not a number, try again. y is','s'));
            end
        end
        
        segment  =   str2double( input('How many segment per edge? Please enter a valid integer','s'));
        while isnan(segment) || fix(segment)~= segment || segment<=1
            segment  =   str2double( input('How many segments per edge? Please enter an positive integer greater than 1','s'));
        end
        
        %%get the friction coefficient
        u = str2double( input(' friction coefficient is','s')    );
        while isnan( u)  || u<0 || u>1
            u = str2double( input('Not a number, try again. friction coefficient is','s')    );
        end
        my_unit = input('what is the unit\n','s');
    end
end