function [M, W1, W2] = CENTERLINE_pombe(x,y, center, stp)
% Complete function for producing centerline and width points 
% of S. pombe contour
%
% 1. Requires input of x,y, and center vectors of contour
% 2. Important aspect is step size, 'stp'
%       - All calculations are in pixels
%       - Recommended step size is ~15 pixels at 100X. Larger
%           seems to work better
%
% Output is M = array of midpoint vectors, W1 = first width
% vectors, W2 = second width vectors
%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Initial midpoint analysis for determining centerline 
%    directions

% Initial width to start from
[w1c, w2c, widc] = calc_wid(x,y,center);

% Obtain vector for width markers
nc = (w2c -w1c)/norm(w2c-w1c);


midc = w1c + nc*0.5*norm(w2c-w1c);


nvc = transpose(null(nc));


for j=1:1
    
    if j==1
        [w1,w2,wid] = calc_wid(x,y,midc+nvc*stp);
        nv = nvc;
        
    else
        
        [w1,w2,wid] = calc_wid(x,y,mid+nv*stp);
        
    end
    
    % Obtain vector for width markers
    n = (w2 -w1)/norm(w2-w1);
    
    mid = w1 + n*0.5*norm(w2-w1);    
        
    
    nvnew = transpose(null(n));
    
    norm(nvnew+nv);
    dot(nvnew,nv);
    if dot(nvnew,nv)<0
        nv = -nvnew;
    else
        nv = nvnew;
    end
    

    
    
    
    
end



nvc1 = nvc;
nvc2 = -nvc;





%%%%%%%%%%%%%%%%%%%
% Store width and midpoint markers

M = []; W1 = []; W2 = [];

M = [M; midc];
W1 = [W1; w1c];
W2 = [W2; w2c];

% Store directional vectors
n1 = nvc1;
n2 = nvc2;




%%%%%%%%%%%%%%%%%%%%%
% 2. Start moving in n1 direction
n = n1;

% Define initial midpoint 
m1 = M(end,:);

% Move to next possible midpoint with vector step
m2 = m1 + n*stp;

% Test if m2 is inside contour
inout = inpolygon(m2(1), m2(2), x,y);

while inout==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Center contour to new midpoint 'mid'

    xn = x - m2(1);
    yn = y - m2(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. Rotate contour from previous midpoint, m1

    m1n = m1 - m2;

    theta = pi/2 -abs(atan(m1n(2)/m1n(1)));

    directions = {'clockwise', 'anticlockwise'};

    if m1n(2)*m1n(1)<0 % tests for being in 2nd or 4th quadrant
        direction = char(directions(1));
    elseif m1n(2)*m1n(1)>0 % tests for 1st, 3rd quadrant
        direction = char(directions(2));
    end

    [xr,yr,xor,yor,temp] = rotateData(xn,yn,0, 0,theta,direction);

    xr = xr';
    yr = yr';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort yr by points greater than zero, and find two smallest values

    % need to sort xr, yr according to yr

    v = [xr yr];

    in = v(:,2)>0;

    vn = v(in,:);

    [ys, is] = sort(vn(:,2));
    vs = vn(is, :);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 6. Find new width markers and set new midpoint

    % Now new width markers are vs(1,:) and vs(2,:)

    w1 = vs(1,:); % Store these in struct
    w2 = vs(2,:);

    dw = norm(w2-w1); % store in struct?

    midn = w1 + (w2-w1)/2; % Defining new midpoint. 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Store in struct as centerline? Will need to rotate back.

    addingstruct = true;

    if addingstruct
        % Rotate back new midpoint
        %if direction==char(directions(1))
        if strcmp(direction, directions(1))
            ndirection = char(directions(2));
        else 
            ndirection = char(directions(1));
        end
        [midx,midy,xor,yor,temp] = rotateData(midn(1),midn(2),0, 0,theta,ndirection);
        [w1x,w1y,xor,yor,temp] = rotateData(w1(1),w1(2),0, 0,theta,ndirection);
        [w2x,w2y,xor,yor,temp] = rotateData(w2(1),w2(2),0, 0,theta,ndirection);
         
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        midcnn = [midx, midy] + m2;
        w1n = [w1x, w1y] + m2;
        w2n = [w2x, w2y] + m2;
        
        
        % This is the vector to place as centerline point 
        
        M = [M; midcnn];
        W1 = [W1; w1n];
        W2 = [W2; w2n];
     
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    
    % Create new vector
    w1 = W1(end,:); w2 = W2(end,:);
    nnew = (w2 -w1)/norm(w2-w1);

    n= transpose(null(nnew));

    % Align vector to initial vector chosen from centroid
    if dot(n, n1) < 0
        n = -n;
    end
    
    
%     % Arrange width vectors
%     dw = norm(w2-w1);
%     dw1 = norm(W1(end,:)-W1(end-1,:));
%     if dw1>dw/2
%         temp1 = W1(end,:);
%         temp2 = W2(end,:);
%         W2(end,:) = temp1;
%         W1(end,:) = temp2;
%     end
    
    
    % Arrange width vectors
    dw = W2(end,:)-W1(end,:);
    dw1 = W2(end-1,:)-W1(end-1,:);
    if dot(dw,dw1)<0
        temp1 = W1(end,:);
        temp2 = W2(end,:);
        W2(end,:) = temp1;
        W1(end,:) = temp2;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%
    % Now rename vectors for looping
    
    % Define initial midpoint 
    m1 = M(end,:);

    % Move to next possible midpoint with vector step
    m2 = m1 + n*stp;
    
    inout = inpolygon(m2(1), m2(2), x,y);
    %inout = 0;
    
end

% Since loop ended with inout = 0, find smallest distance point
% to contour
for k=1:length(x)
    dd(k, :) = norm(m2 - [x(k) y(k)]);
end

[ddmin, I] = min(dd);

mend1 = [x(I) y(I)];

M = [M; mend1];

% Also assign midpoint as width markers as well
w1 = mend1; w2 = mend1;

W1 = [W1; w1];
W2 = [W2; w2];





%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
% 3. Start moving in n2 direction
n = n2;

% Define initial midpoint 
m1 = M(1,:);

% Move to next possible midpoint with vector step
m2 = m1 + n*stp;

% Test if m2 is inside contour
inout = inpolygon(m2(1), m2(2), x,y);

% Make new storing vectors

M_2 = M(1,:);
W1_2 = W1(1,:);
W2_2 = W2(1,:);



while inout==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Center contour to new midpoint 'mid'

    xn = x - m2(1);
    yn = y - m2(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. Rotate contour from previous midpoint, m1

    m1n = m1 - m2;

    theta = pi/2 -abs(atan(m1n(2)/m1n(1)));

    directions = {'clockwise', 'anticlockwise'};

    if m1n(2)*m1n(1)<0 % tests for being in 2nd or 4th quadrant
        direction = char(directions(1));
    elseif m1n(2)*m1n(1)>0 % tests for 1st, 3rd quadrant
        direction = char(directions(2));
    end

    [xr,yr,xor,yor,temp] = rotateData(xn,yn,0, 0,theta,direction);

    xr = xr';
    yr = yr';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort yr by points greater than zero, and find two smallest values

    % need to sort xr, yr according to yr

    v = [xr yr];

    in = v(:,2)>0;

    vn = v(in,:);

    [ys, is] = sort(vn(:,2));
    vs = vn(is, :);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 6. Find new width markers and set new midpoint

    % Now new width markers are vs(1,:) and vs(2,:)

    w1 = vs(1,:); % Store these in struct
    w2 = vs(2,:);

    dw = norm(w2-w1); % store in struct?

    midn = w1 + (w2-w1)/2; % Defining new midpoint. 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Store in struct as centerline? Will need to rotate back.

    addingstruct = true;

    if addingstruct
        % Rotate back new midpoint
        %if direction==char(directions(1))
        if strcmp(direction, directions(1))
            ndirection = char(directions(2));
        else 
            ndirection = char(directions(1));
        end
        [midx,midy,xor,yor,temp] = rotateData(midn(1),midn(2),0, 0,theta,ndirection);
        [w1x,w1y,xor,yor,temp] = rotateData(w1(1),w1(2),0, 0,theta,ndirection);
        [w2x,w2y,xor,yor,temp] = rotateData(w2(1),w2(2),0, 0,theta,ndirection);
         
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        midcnn = [midx, midy] + m2;
        w1n = [w1x, w1y] + m2;
        w2n = [w2x, w2y] + m2;
        
        
        % This is the vector to place as centerline point 
        
        M_2 = [M_2; midcnn];
        W1_2 = [W1_2; w1n];
        W2_2 = [W2_2; w2n];
     
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    
    % Create new vector
    w1 = W1_2(end,:); w2 = W2_2(end,:);
    nnew = (w2 -w1)/norm(w2-w1);

    n= transpose(null(nnew));

    % Align vector to initial vector chosen from centroid
    if dot(n, n2) < 0
        n = -n;   
    end
    
    
    % Arrange width vectors
    dw = W2_2(end,:)-W1_2(end,:);
    dw1 = W2_2(end-1,:)-W1_2(end-1,:);
    if dot(dw,dw1)<0
        temp1 = W1_2(end,:);
        temp2 = W2_2(end,:);
        W2_2(end,:) = temp1;
        W1_2(end,:) = temp2;
    end
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%
    % Now rename vectors for looping
    
    % Define initial midpoint 
    m1 = M_2(end,:);

    % Move to next possible midpoint with vector step
    m2 = m1 + n*stp;
    
    inout = inpolygon(m2(1), m2(2), x,y);
    %inout = 0;
    
end

% Since loop ended with inout = 0, find smallest distance point
% to contour
for k=1:length(x)
    dd(k, :) = norm(m2 - [x(k) y(k)]);
end

[ddmin, I] = min(dd);

mend1 = [x(I) y(I)];

M_2 = [M_2; mend1];

% Also assign midpoint as width markers as well
w1 = mend1; w2 = mend1;

W1_2 = [W1_2; w1];
W2_2 = [W2_2; w2];



%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
% Now invert M_2, W1_2, W2_2 and catenate with original storing
% vectors!

M = [flipud(M_2(2:end,:)); M];
W1 = [flipud(W1_2(2:end,:)); W1];
W2 = [flipud(W2_2(2:end,:)); W2];

% M = [M; flipud(M_2(2:end,:))];
% W1 = [W1; flipud(W1_2(2:end,:))];
% W2 = [W2; flipud(W2_2(2:end,:))];













