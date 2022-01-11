%1. Random 3x4 matrix
    %The for loops create a 3x4 matrix where each element is randomly
    %generated
    for col = 1:4
        for row = 1:3
            rmat(row,col) = rand;
        end    
    end
    disp(rmat);
    %Identifies the smallest value in each row where M is a column vector
    %with each of the smallest value and I is the column that the smallest
    %value in each row appears in
    [M,I] = min( rmat, [], 2 );
    disp('Smallest values per row')
    disp(M)
    disp ('Column that smallest value appears in')
    disp(I)

%Displays car information
    car = struct('color',' white','year',' 1998','make',' Honda','model',' Civic','class',' LX','cylinders',' 4 cylinders');
    disp(['Say you drive a ',  car.color, car.year , car.make, car.model ,  car.class , car.cylinders ])
    
%Displays car information but uses a cell
    carcell = {'color','year','make','model','class','cylinder';' white',' 1998',' Honda',' Civic',' LX',' 4 cylinders'};
    disp(['Say you drive a ', carcell(2,1), carcell(2,2), carcell(2,3), carcell(2,4), carcell(2,5), carcell(2,6)])

%Make a 10Hz sine wave
    t = 0:0.2:0.6;
    y = 3*sin(2*sin(2*pi*10*t));
    plot(t,y)
    %It doesn't look like a sine wave because there are to few steps between
    %the starting and ending t values

%Make another sine wave
    t = 0:0.001:100; 
    y = 3*sin(2*sin(2*pi*10*t));
    plot(t,y)
    %The window is too zoomed and the lines of the sine function all run
    %together

%An actual sine wave
    %The first and last values of f represent the beginnning and the end while
    %the second is the size of sections of t and determines the resolution of
    %the graph
    t = 0:0.001:.2; 
    y = 3*sin(2*sin(2*pi*10*t));
    plot(t,y,'LineWidth',4);
   %Formatting the plot
    grid
    xlabel('time');
    ylabel('height');
    title('Sine Wave');
   


