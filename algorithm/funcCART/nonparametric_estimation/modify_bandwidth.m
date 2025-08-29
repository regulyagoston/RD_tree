%%%
% Modifies the bandwidth value dependent on the sample size
% Following Calanico et al (2020, EJ, p200)

function [ bw ] = modify_bandwidth( bwselect, orig_bw, orig_n, new_n, ord_poly )

switch bwselect
    case {"mserd", "msesum","msecomb1"}
        exp_val = -1/(3+2*ord_poly);
    case {"cerrd" ,"cersum","cercomb1"}
        exp_val = -1/(3+ord_poly);
    otherwise
        error('Not implemented')
end

bw = orig_bw .* (new_n./orig_n).^exp_val;

end