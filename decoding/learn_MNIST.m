images = loadMNISTImages('train-images-idx3-ubyte'); % initialize figure  
labels = loadMNISTLabels('train-labels-idx1-ubyte'); % initialize figure
labels = labels';                                    % transpose
labels(labels==0)=10;                                % dummyvar function doesn´t take zeroes
labels=dummyvar(labels);  %

% get a feeling of the images
get(groot,'default')
figure                                          % initialize figure
colormap(gray)                                  % set to grayscale
for i = 1:36                                    % preview first 36 samples
    subplot(6,6,i)                              % plot them in 6 x 6 grid
    digit = reshape(images(:, i), [28,28]);     % row = 28 x 28 image
    imagesc(digit)                              % show the image
    title(num2str(labels(i)))                   % show the label
end


% construction of neural net
x = images;
t = labels';
trainFcn = 'trainscg';                          % use scaled conjugate gradient for training

hiddenLayerSize = 100;                          
net = patternnet(hiddenLayerSize);              % create Pattern Recognition Network

net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
net.performFcn = 'crossentropy';

[net,tr] = train(net,x,t);
