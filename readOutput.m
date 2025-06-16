Output_param = readmatrix("ModelOutput.csv");

smallest_optimised_param_costfunction = min(Output_param(15,:));
disp(smallest_optimised_param_costfunction)
col_num_of_smallest_op_param = find(Output_param(15,:)==smallest_optimised_param_costfunction);
optimised_param = Output_param(:, col_num_of_smallest_op_param).'; %this matrix would include the cost function at the end for reference
disp(optimised_param)