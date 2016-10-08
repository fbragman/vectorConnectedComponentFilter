classdef EquivalencyTable < handle
    %  EquivalencyTable is based on itkEquivalencyTable. It is a hash table
    %  for recording equivalencies among unsigned integers.
    %  EquivalencyTable can store recursive relationships (8=7, 7=6, 6=5)
    %  or be flattened to eliminate recursion (8=5,7=5,6=5)
    %
    %  Written by Felix Bragman
    %  UCL Centre for Medical Image Computing (2015). 
    
    % Private class properties
    properties (SetAccess = private)
        %% Private properties
        equivalencyTable;
        flatTable;
    end
          
    methods
        
        function obj = EquivalencyTable()
           %%   Equivalency table constructor 
            
           % Constructor. Create an equivalency table based on the number
           % of voxels to analyse and the dimensionality of the image           
           eqTable = [];
           obj.equivalencyTable = eqTable;
        end
        
        
        
        function Add(obj,a,b)
            %%  Add equivalency to table
            
            % Insert an equivalency into the table.
            
            table = obj.equivalencyTable;
            tmp   = [a b];
            
            if isempty(table)
                % if nothing has been added yet
                table = tmp;
            else
                % append the table then sort
                table(end+1,:) = tmp;
                [~,idxSort] = sort(table(:,1),'ascend');
                table = table(idxSort,:);
            end
            
            obj.equivalencyTable = table; 
        end
        
        function result = VectorisedLookUp(obj,vec)
           %%   Vectorised Lookup 
           
           % Extract object
           fTable = obj.flatTable;
           
           % Equivalency of labels
           [idxA,locA] = ismember(vec,fTable(:,1));
           locA(locA==0) = [];
           
           result = zeros(size(vec));
           % If not found in flattened equivalency table, then the label
           % remains the same
           result(~idxA) = vec(~idxA);
           % If found in flattened equivalency table, set label to
           % flattened label
           result(idxA)  = fTable(locA,2);
        end
        
        function result = LookUpRecursive(obj,start)
            %%  Recursive Look Up
            
            table = obj.equivalencyTable;
            
            if isempty(table)
                
                result = start;
                
            else
            
                % Initialise recursive look-up with starting value a
                lookedAtNodes = start;
                
                % automatic growth of lookedAtNodes for optimisation
                
                % Initialise while loop exit flag
                flag = 1;
                
                % starting from a parent node e.g 10, look at all equivalences going
                % down/up the connected graph, finding the lowest connection at each
                % point and stopping when no more connections are found
                while (flag)
                    
                    % end loop prematurely if there is only one occurence of the value
                    % allowing one to immediately go to the next node in the graph.
                    %
                    % we do not operate on minValues to guard against cases where we might
                    % miss more complicated node connections
                    
                    existInMin = ismembc(table(:,1),start);
                    if ~existInMin
                        break
                    end
                    
                    [tmp,numU] = count_unique(table(existInMin,1));
                    if ~any( ones(numel(tmp),1) - numU )
                        start = table(existInMin,2);
                        lookedAtNodes = [lookedAtNodes start'];
                        continue
                    end
                    
                    attached2Parent = table(ismembc(table(:,1),start),2);
                    
                    if isempty(attached2Parent)
                        % no equivalency (yet)
                        break
                    end
                    
                    % guard against going back up the tree
                    % e.g. [ 7 ..
                    %        8 10
                    %        10 7
                    %        10 8 ]
                    attached2Parent = uint32(count_unique(attached2Parent));
                    checkDone = ismembc(attached2Parent,lookedAtNodes);
                    if all(checkDone)
                        break
                    end
                    
                    % add successive nodes into lookedAtNodes, which do not complete loops
                    % e.g. [ 8  10
                    %        10 8 ]
                    %
                    lookedAtNodes = [lookedAtNodes attached2Parent(~checkDone)'];
                    if numel(find(ismembc(attached2Parent,table(:,1)))) == 0
                        break
                    else
                        start = attached2Parent(~checkDone);
                    end
                end
                result = min(lookedAtNodes);
            end
        end
        
        function result = LookUp(obj,a)
           %%   Table look up
            
           % Look up equivalency in table. If no entry is found in the
           % table, the method will return the value of the argument. Not a
           % recursive lookup like in LookUpRecursive. Operated on a
           % flattened table
           
           fTable = obj.flatTable;
           
           if isempty(fTable)
               result = a;
           else
               [locA,locB] = ismember(a,fTable(:,1));
               
               if ~locA
                   result = a;
               else
                   result = fTable(locB,2);
               end
           end
        end
        
        function Flatten(obj)
           %%   Flatten table
           % Flatten equivalency table using the construction of an
           % adjancency matrix and Dulmage-Mendelsohn permutation
           
           table = obj.equivalencyTable;
           % Create proper pairing 1 <=> 2 --> [1 2;2 1]
           tmp2 = table;
           tmp3 = table;
           tmp2(:,1) = tmp3(:,2);
           tmp2(:,2) = tmp3(:,1);
           clear tmp3;
           masterT = zeros(2*size(table,1),2);
           masterT(1:2:end-1,:) = table;
           masterT(2:2:end,:)   = tmp2;
           clear tmp2
           masterT = double(masterT);
           
           % Create of adjacency matrix
           A = masterT;
           spA = sparse(A(:,1),A(:,2),1);
           spA = spA + speye(size(spA,1),size(spA,2));
           
           % Dulmage-Mendelsohn permutation
           [p,~,r,~] = dmperm(spA);
           
           % Extracting of equivalent labelling components
           cellRes = cell(numel(r)-1,1);
           for idx = 1:numel(r)-1
               cellRes{idx} = p(r(idx):r(idx+1)-1);
           end
           
           % Deletion of single components (?)
           id = cellfun(@(x) numel(x) == 1,cellRes);
           cellRes(id) = [];
           
           % Creation of component matrices
           tmp = cell(size(cellRes));
           for idx = 1:numel(cellRes)
               tmpM    = cellRes{idx};
               tmpMin  = min(tmpM);
               tmpM(tmpM==tmpMin) = [];
               mat = [tmpM(:) repmat(tmpMin,numel(tmpM),1)];
               tmp{idx} = mat;
           end
           
           % Conversion to table
           flat = cell2mat(tmp);
           
           obj.flatTable = uint32(flat);
        end
            
        function result = Empty(obj)
           %%   Table empty check 
            
           % Returns TRUE if the table is empty. FALSE if it is not empty. 
            
           if isempty(obj.equivalencyTable) 
               result = 1;    
           else    
               result = 0; 
           end  
        end
        
        function SetEqTable(obj,table)
           %%   Setter
           
           obj.equivalencyTable = table;
            
        end
        
    end     
end  

function result = innerLookUp(table,lookup)
%%  Inner look up function for Flatten

[locA,locB] = ismember(lookup,table(:,1));

if ~locA
    result = lookup;
else
    result = table(locB,2);
end

end