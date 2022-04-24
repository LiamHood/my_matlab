close all; clear; clc;
warning('off','all')
iris_tests([4 3], 'traingd', 'logsig')
iris_tests([2 3], 'traingd', 'logsig')
iris_tests([1 3], 'traingd', 'logsig')
iris_tests([4 3], 'traingd', 'tansig')
iris_tests([4 3], 'trainlm', 'logsig')
iris_tests([4 3], 'trainlm', 'tansig')
function iris_tests(architecture, learn, act)
    fprintf("Hidden: %f\tLearning Rule: %s\tActivation Function: %s\n", ...
        [architecture(1), learn, act])
    load("iris.mat")
    iris_tra_input = iris_tra(:, 1:4)';
    iris_tra_target = iris_tra(:, 5:7)';
    iris_tes_input = iris_tes(:, 1:4)';
    iris_tes_target = iris_tes(:, 5:7)';
    
    netiris1 = newff(minmax(iris_tra_input), architecture, {act, act}, learn);
    
    netiris1.trainParam.show = 50;
    netiris1.trainParam.lr = .05;
    netiris1.trainParam.epochs = 10000;
    netiris1.trainParam.goal = 1e-5;
    [netiris1, tr] = train(netiris1, iris_tra_input, iris_tra_target);
    perf_tra_1 = perform(netiris1, iris_tra_target, netiris1(iris_tra_input));
    perf_tes_1 = perform(netiris1, iris_tes_target, netiris1(iris_tes_input));
    rms_tra_1 = sqrt(mean(mean((netiris1(iris_tra_input) - iris_tra_target).^2)));
    rms_tes_1 = sqrt(mean(mean((netiris1(iris_tes_input) - iris_tes_target).^2)));
    
    netiris2 = feedforwardnet(architecture(1), learn);
    netiris2 = configure(netiris2, iris_tra_input, iris_tra_target);
    netiris2.layers{1}.transferFcn = act;
    netiris2.layers{2}.transferFcn = act;
    netiris2.trainParam.show = 50;
    netiris2.trainParam.lr = .05;
    netiris2.trainParam.epochs = 10000;
    netiris2.trainParam.goal = 1e-5;
    [netiris2, tr] = train(netiris2, iris_tra_input, iris_tra_target);
    bb = netiris2(iris_tra_input);
    perf_tra_2 = perform(netiris2, iris_tra_target, netiris2(iris_tra_input));
    perf_tes_2 = perform(netiris2, iris_tes_target, netiris2(iris_tes_input));
    rms_tra_2 = sqrt(mean(mean((bb - iris_tra_target).^2)));
    rms_tes_2 = sqrt(mean(mean((netiris2(iris_tes_input) - iris_tes_target).^2)));
    
    fprintf("newff\nRMS Training: %f\tRMS Testing: %f\tPerf Training: %f\tPerf Testing: %f\n", ...
        [rms_tra_1, rms_tes_1, perf_tra_1, perf_tes_1])
    fprintf("feedforward\nRMS Training: %f\tRMS Testing: %f\tPerf Training: %f\tPerf Testing: %f\n\n", ...
        [rms_tra_2, rms_tes_2, perf_tra_2, perf_tes_2])
end