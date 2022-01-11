clear ; close all ; clc ;
% intitial set up
for ii = 1:3
    for jj = 1:3
        squares{ii,jj} = zeros(3,3) ;
    end
end
board = [ squares{1,1} , squares{1,2} , squares{1,3} ; squares{2,1} , squares{2,2} , squares{2,3} ; squares{3,1} , squares{3,2} , squares{3,3} ] ;

% assigning numbers
kk = 1 ;
while 0 == isempty( find( board == 0 ) )
            
            % finds first 0 on the board to replace with a number
            [ ii , jj ] = find( board == 0 ) ;
            ii = ii( 1 ) ;
            jj = jj( 1 ) ;
            % picks a random number between 1 and 9
            num = randi(9) ;
                    % checks for row repeat
                    checkr1 = board( ii , 1:9 ) ;
                    checkr = find( num == checkr1 ) ;
                    % checks for column repeat
                    checkc1 = board( 1:9 , jj ) ;
                    checkc = find( num == checkc1 ) ;
                    % checks for repeat in square
                    checkb1 = squares{ ceil(ii/3),ceil(jj/3) }  ;
                    checkb = find( num == checkb1 ) ;
                    if isempty( checkr )
                        if isempty( checkc )
                            if isempty( checkb )
                                % if no repeats the number is assigned to
                                % the cell and resets attempt counter
                                board(ii,jj) = num ; 
                                kk = 1 ;
                            else 
                                % failure of any repeat results in an
                                % addition to the attempt counter
                                kk = kk + 1 ;
                            end
                        else 
                            kk = kk + 1 ;
                        end
                    else 
                            kk = kk + 1 ;    
                    end
            % adds the new number to both the board and the squares
            squares{ceil(ii/3),ceil(jj/3)} = board( (ceil(ii/3)*3-2):(ceil(ii/3)*3) , (ceil(jj/3)*3-2):(ceil(jj/3)*3) ) ;
            board = [ squares{1,1} , squares{1,2} , squares{1,3} ; squares{2,1} , squares{2,2} , squares{2,3} ; squares{3,1} , squares{3,2} , squares{3,3} ] ;
    % Reset board if a number can't be found that checks out. kk counts times that a random number has been assigned     
    if kk == 20
        for ii = 1:3
            for jj = 1:3
                squares{ii,jj} = zeros(3,3) ;
            end
        end
        board = [ squares{1,1} , squares{1,2} , squares{1,3} ; squares{2,1} , squares{2,2} , squares{2,3} ; squares{3,1} , squares{3,2} , squares{3,3} ] ;
    end
end
