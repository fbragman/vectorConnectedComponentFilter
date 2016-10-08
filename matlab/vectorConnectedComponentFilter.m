function [result,eqTable] = vectorConnectedComponentFilter(img,mask)
%%  Vector Connected Component Image Filter
%
%   Description: vectorConnectedComponentFilter is an implementation of the
%   itkVectorConnectedComponentImageFilter. It labels components in an
%   image based on a pre-determined face-connectivity (4/6 2D/3D) and a
%   function (dot product of 2 vectors). vectorConnectedComponentFilter
%   returns a labeled image of the components, which can be used for
%   further processing. Algorithm hard-coded for 6 connectivity 3D images
%   vector images.
%
%   Syntax: result = vectorConnectedComponentFilter(img,mask)
%
%   Input:  img     -   vector image
%           mask    -   mask of vector image
%
%   Output: result  -   labelled image of img
%
% Written by Felix Bragman, UCL Centre for Medical Image Computing (2015).

fprintf('Running vcc \n');

% Similarity value
s = 0.90;
sim = 1-s;

% Dot product functor
simFoo = @(x,y,S) 1 - abs(x*y') <= S;

% Dimension
d = 3;

% Algorithm
[result,eqTable] = VCCfoo(img,mask,d,sim,simFoo);
fprintf('Finished vcc \n');
end

function [output,eqTable] = VCCfoo(img,mask,d,S,functor)
%%  Test algorithm
%
%   Inputs:     img     -   vector image
%               conn    -   connectivity
%               d       -   image dimension
%
%   Operate on padded mask with padded value given label 0

padVector = ones(1,d);

% Img vector
vec1     = img(:,:,:,1);
vec1_pad = padarray(vec1,padVector,0,'pre'); clear vec1; 
vec2     = img(:,:,:,2);
vec2_pad = padarray(vec2,padVector,0,'pre'); clear vec2;
vec3     = img(:,:,:,3);
vec3_pad = padarray(vec3,padVector,0,'pre'); clear vec3;

clear img;

% Pad binary mask
pad_mask = padarray(mask,padVector,0,'pre');
siz_pad  = size(pad_mask);

% Changed - hard-code x+1,y+1,z+1 instead
% Linear index conversion 
% oit_conv = convertIterator(siz_pad);

% Index image
[m,n,b]   = size(pad_mask);
N = m*n*b;

% MAIN ENGINE: iterate over the image, labelling the objects and defining
% equivalence. Use the connectivity (neighbourhood iteration) to access
% previously analysed neighbour pixels and an output iterator to access the
% current pixel

% initialise label to maximum allowable by the data class
labelType = 'uint32';
maxPossibleLabel  = intmax(labelType);
label_img = reshape(maxPossibleLabel*ones(1,N,labelType),m,n,b);

% use mask to mark pixel as unlabelled (same as the padded voxels to be
% ignored)
label_img(~pad_mask) = 0;
% label_img = padarray(label_img,padVector,0,'pre');

% initialise label count
maxLabel = 0;

% initialise equivalency table
% column 1 - current node being traversed
% column 2 - neighbour 1 (0 -1)
% column 3 - neighbour 2 (-1 0)
eqTable = EquivalencyTable();

% C++ this loop
fprintf('Starting loop \n');
        
    idxInImg    = find(label_img);
    % Subscript (non-padded)
    mask_subV   = ind2subv(size(label_img),idxInImg);
    % Sort such that we start at (1,1,1) and go down in following fashion
    %
    % | o o x x |     | 0 0 11 16 |     
    % | x o x x |     | 2 0 12 17 |
    % | x o x x | --> | 3 0 13 18 | --> [2,3,9,10,11,12,13,16,17,18,19,20]
    % | o x o x |     | 0 9  0 19 |
    % | o x o x |     | 0 10 0 20 |
    mask_subVso = sortrows(mask_subV,[3 2 1]);
    idxInImgso  = subv2ind(size(label_img),mask_subVso);
tic    
    % idxInImg  = find(label_img);
    for idx = 1:numel(idxInImgso)
        fprintf('%f percent complete \n',(idx/numel(idxInImgso))*100);
 
        it  = idxInImgso(idx);
        oit = idxInImgso(idx);
        
        % Directly converted
%         label = label_img(oit_conv(it));
        label = label_img(it);
        
%         value = [vec1_pad(oit_conv(it)) ...
%             vec2_pad(oit_conv(it)) ...
%             vec3_pad(oit_conv(it))];
        
        % operating on logical index of padded array, not need to convert
        value = [ vec1_pad(it) ...
                  vec2_pad(it) ...
                  vec3_pad(it) ];
        
        originalLabel = label;
        
        % Find previous connected neighbour
%         [init,onit] = nhoodIterator(oit_conv(oit),siz_pad);
        [init,onit] = nhoodIterator(oit,siz_pad);
        
        for niter = 1:d
            neighbourLabel = label_img(onit(niter));
            
            % Check that neighbourhood voxel is not a padded voxel or
            % background
            if (logical(neighbourLabel))
                
                neighbourValue = [vec1_pad(init(niter)) ...
                    vec2_pad(init(niter)) ...
                    vec3_pad(init(niter))];
                
                % Dot product of normalised eigenvectors
                if (functor(value,neighbourValue,S))
                    
                    % If current voxel is unlabelled, copy label from neighbour
                    if (label == maxPossibleLabel)
                        
                        % Copy label from previous pixel
                        label = neighbourLabel;
                        
                        % Else if current pixel has a label that is not already
                        % equivalent to the label of the previous pixel, then setup
                        % a new equivalence
                        
                        % Inequality condition:
                        % If label doesn't equal neighbour label and we know they
                        % aren't yet connected by recursively looking at
                        % connections of both labels
                    elseif label ~= neighbourLabel && ...
                            ( eqTable.LookUpRecursive(label) ~= ...
                            eqTable.LookUpRecursive(neighbourLabel) )
                        
                        if neighbourLabel > label
                            
                            eqTable.Add(neighbourLabel,label);
                            
                        else
                            
                            eqTable.Add(label,neighbourLabel);
                            
                        end
                        
                    end
                end
            end
        end
        
        % If none of the "previous" neighbours were set, make a new label
        if (originalLabel == label)
            
            % Create a new entry label
            maxLabel = maxLabel + 1;
            label = maxLabel;
        end
        
        % Set output pixel to label we currently have
        if (label ~= originalLabel)
            
%             label_img(oit_conv(oit)) = label;
            label_img(oit) = label;
            
        end
        
        % oit   = oit + 1;
        % it    = it + 1;
    end
toc

% Flatten equivalency table
eqTable.Flatten();
% % Equivalence class resolution second pass
idx_label     = find(label_img);
tmp           = label_img(idx_label);
newLabels     = eqTable.VectorisedLookUp(tmp);

output = zeros(size(label_img));
output(idx_label) = newLabels;

% Remove label_img padding
output = output(1+padVector(1):end,1+padVector(2):end,1+padVector(3):end);
end

function oit_conv = convertIterator(siz_pad)
%%  Index converter between padded and non-padded image

N_pad = prod(siz_pad);

oit_conv = 1:N_pad;

ind_pad = sub2ind(siz_pad,2,2);
tmp_ind = cumsum([1 repmat(siz_pad(1),1,siz_pad(1)-1)]);

oit_conv(tmp_ind) = [];
oit_conv(oit_conv<ind_pad)  = [];
end

function [init,onit] = nhoodIterator(oit,siz)
%%  Neighbouhood iterator

funSub = @(pos) repmat(pos,3,1) - eye(3);

node_sub = ind2subv(siz,oit);
sub_previous = funSub(node_sub);

previous_idx = subv2ind(siz,sub_previous);

init = previous_idx;
onit = previous_idx;
end