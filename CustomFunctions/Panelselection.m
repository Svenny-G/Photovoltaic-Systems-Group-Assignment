function [selected_modules,num_panels_req,expected_yield] = Panelselection(annual_yield,Req_generation_X)
    [sorted_yield, sorted_indices] = sort(annual_yield, 'descend');

    cumsumyield = cumsum(sorted_yield);
    num_panels_req = find(cumsumyield >= Req_generation_X, 1, "first");
    selected_modules = sorted_indices(1:num_panels_req);
    expected_yield = cumsumyield(num_panels_req);