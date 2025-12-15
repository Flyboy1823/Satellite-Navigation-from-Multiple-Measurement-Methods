function result = EL_AZ_Ploter(el_matrix,az_matrix)
%--------------------------------------------------------------------------
% Generates a sky plot for satellites whose relative elevation and azimuth
% values are given as matrices.
%--------------------------------------------------------------------------
%
% input:    'svid'          a 1 x num_sv or num_sv x 1 vector of SV prn numbers
%
%       	'el_matrix'     a N x 1 array which contains a series of elevation measurements [deg]
%		
%       	'az_matrix'     a N x 1 array which contains a series of azimuth measurements [deg]
%
%
% output: 'result' indicates if the data have been plotted successfully 

colormap(lines);
cmap=colormap;
color_dark_gray=[1,1,1]*0.5;
color_light_gray=[1,1,1]*0.95;
s_vec=linspace(0,2*pi,101);

% Figure construction
figure;

% Create auxiliary axes and marking
plot([-90,90],[0,0],'Color',color_dark_gray);
hold on; axis equal;
plot([0,0],[-90,90],'Color',color_dark_gray);
plot(30*sin(s_vec),30*cos(s_vec),'Color',color_dark_gray);
plot(60*sin(s_vec),60*cos(s_vec),'Color',color_dark_gray);
plot(90*sin(s_vec),90*cos(s_vec),'Color',color_dark_gray);
text(-3,95,'N','Color',color_dark_gray);
text(3,63,'30','Color',color_dark_gray);
text(3,33,'60','Color',color_dark_gray);
text(3,3,'90','Color',color_dark_gray);
text(-3,-95,'S','Color',color_dark_gray);

    
% Check visibility changes in time for current SV
sign_vec=[0;sign(el_matrix(:))];
diff_sign_vec=diff(sign_vec);
diff_sign_vec_ind=[find(diff_sign_vec~=0);size(el_matrix,2)];

% draw visible/hidden trajectory segments 
for index2=1:length(diff_sign_vec_ind)-1
    
    % Calc values for current segment
    sign_type=diff_sign_vec(1)*(-1)^(index2-1);
    el_values=el_matrix(diff_sign_vec_ind(index2):diff_sign_vec_ind(index2+1));
    az_values=az_matrix(diff_sign_vec_ind(index2):diff_sign_vec_ind(index2+1));
    x_values=(90-abs(el_values)).*(sin(az_values*pi/180));
    y_values=(90-abs(el_values)).*(cos(az_values*pi/180));
    
    if sign_type==1
        % Case of a visible segment
        h_line=plot(x_values,y_values,'Linewidth',2);
    else
        % Case of an hidden segment
        h_line=plot(x_values,y_values,'LineStyle','--');
    end
    
    % Additional marking for trajectory start/end points
    set(h_line,'Color',cmap(1,:));
    if index2==1
        h_init=plot(x_values(1),y_values(1),'o','Color',cmap(1,:),'MarkerSize',8);
        %text(x_values(1)+5,y_values(1),mat2str(svid(index1)),'Color',cmap(index1,:),'BackgroundColor',color_light_gray);
        if sign_type==1
            set(h_init,'LineWidth',2);
        end
    end
    if index2==(length(diff_sign_vec_ind)-1)
        h_end=plot(x_values(end),y_values(end),'x','Color',cmap(1,:),'MarkerSize',8);
        if sign_type==1
            set(h_end,'LineWidth',2);
        end
    end
    
end
set(gca,'Xlim',[-100,100],'Ylim',[-100,100],'Xtick',[],'Ytick',[]);
result=1;
end







