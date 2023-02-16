clear 

%%
g = 9.81; % gia toc trong truong
length_x = 200; % do dai cua mien gia tri x
length_y = 200; % do dai cua mien gia tri y
num_x = 201; % so khoang chia cua mien gia tri x
num_y = 201; % so khoang chia cua mien gia tri cua y
division_x = length_x / (num_x - 1); % do chia cua mien x
division_y = length_y / (num_y - 1); % do chia cua mien y
height_wall = 2; % chieu cao cua tuong
r = 0.00001; % toc do tat dan
[x, y] = meshgrid(linspace(0,length_x, num_x), linspace(0, length_y, num_y));
%%
U_t = zeros(num_x, num_y);
U_t_plus = zeros(num_x, num_y);
U_t_sub = zeros(num_x, num_y);
wall = zeros(num_x, num_y);
%% read data from txt file
anycase = input('Nhap truong hop can mo phong = ');
switch anycase
    case 1
        f = fopen('n_nguon_khongtatdan.txt','r');
    case 2
        f = fopen('n_nguon_tatdan.txt','r');
    case 3
        f = fopen('n_nguon_tatdan_covatcan.txt','r');
    case 4
        f = fopen('muaroi.txt','r');
    otherwise
        fprintf('\nFile input khong ton tai!\nVui long chon lai!\n');
end
wave_source_num = fscanf(f,'SOURCE_NUM %d\n');
wall_num = fscanf(f,'WALL_NUM %d\n');
damped = fscanf(f,'DAMPED %d\n');
rain = fscanf(f,'RAIN %d\n');
if(rain)
    step = fscanf(f, '%d\n');
end

% doc toa do cua cac nguon song da cho
k = 2;
for num = 1:wave_source_num
    temp = fscanf(f,'WS_ID %d %d %d\n',[3 1]);
    
    % obj_id là tọa độ các nguồn sóng trong mô phỏng
    % với x0 là hoành độ và y0 là tung độ của nguồn sóng được khởi tạo ban 
    x0 = temp(2);
    y0 = temp(3);
    for i=1:num_x 
        for j=1:num_y
            if (wall(i,j))
                U_t_sub(i,j) = 0;          % tường thì để là 0
            else
                U_t_sub(i,j) = U_t_sub(i,j) + 10*exp((-((i-x0)^2 + (j-y0)^2))/(k^2)); 
                
            end
       
        end
    end
end

for num = 1:wall_num
    temp = fscanf(f,'WALL_ID %d %d %d %d %d\n',[5 1]);
    % đọc thông tin về vị trí của các bức tường trong mô phỏng
    type = temp(2);
    u = temp(3);
    v = temp(4);
    w = temp(5);
    if (type == 0)
        wall(v:w,u) = height_wall; % height độ cao của tường
    else
        wall(u,v:w) = height_wall;
    end
end
U_t = U_t_sub;

surf(x, y, U_t)
view(30, 90)
delta_t = 1/30;
Num_steps = 5000; % so buoc lap cua qua trinh mo phong
rt = 0; % he so tat dan
damp = 1;
accT = 0;
n = 0;
speed = 3; % toc do truyen song
ds2 = 0.1;
while n < Num_steps
    % neu co tat dan thi thuc hien theo cong A = A_0 * e^-rt
    if(damped)
        if(rt<2e-2)
            rt = rt + delta_t * r; 
            damp = exp(-rt);
        end
    end
    % ????
    for i = 2:num_x - 1
        if(not(wall(i,1)))
            U_t(i,1) = U_t(i,3);
        end
        if(not(wall(i,num_y)))
            U_t(i,num_y) = U_t(i,num_y - 2);
        end
    end
    for i = 2:num_y-1
        if(not(wall(1,i)))
            U_t(1,i) = U_t(3,i);
        end
        if (not(wall(num_x,i)))
            U_t(num_x,i) = U_t(num_x-2,i);
        end
    end
    
    % tinh toan dao ham bac hai bang phuong phap so
    for i = 2:num_x-1
        for j = 2:num_y-1
            if (wall(i,j))
                continue;
            end
            U_t_plus(i,j)=(2*U_t(i,j)-U_t_sub(i,j)+(ds2)*...
            (U_t(i+1,j)+U_t(i-1,j)+U_t(i,j+1)+U_t(i,j-1)-4*U_t(i,j)))*damp;
			% damp là biến tắt dần
			% công thức trên là ct trong slide 
			% etap là 𝑢(t + ∆t) , eta là u(t) , etam là 𝑢(t - ∆t)
        end
    end
    U_t_sub = U_t;
    U_t = U_t_plus;
    
    accT = accT + 1;
    if(accT < speed)
       continue;
    end
    accT = 0;
    % truong hop neu chon mo phong che do mua roi
    if (rain && mod(n,step)==0)  % chế độ mưa rơi
        while 1
            x0 = randi([1 num_x],1,1);    % thêm các nguồn mới ngẫu nhiên
            y0 = randi([1 num_y],1,1);  
            if (not(wall(x0,y0))) 
                break;
            end
        end
        for i = 1:num_x 
            for j = 1:num_y
                if (wall(i,j))
                    U_t(i,j) = 0;
                else
                    dist = (i-x0)^2+(j-x0)^2;                
                    U_t_sub(i,j) = U_t_sub(i,j) - 5*exp((-((i-x0)^2 + (j-y0)^2))/(k^2));
                    U_t(i,j) = U_t(i,j) - 5*exp((-((i-x0)^2 + (j-y0)^2))/(k^2));
                end
            end
        end
    end

    f = figure(1);
    u = U_t + wall;
    surf(x,y,u)
    colormap(parula(5))
    zlim([-5 5])
    view(45 ,60)
    shading(gca,'interp')
    n = n + 1;
end

        
