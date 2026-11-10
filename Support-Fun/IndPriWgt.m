function w_ip = IndPriWgt(y,Lth)
Nite = length(y);
w_ip = cell(Nite,1);
for ite = 1:Nite
    Nsze = size(y{ite});
    Nsze = Nsze(2:end);
    ind = y{ite}(3,:)>=Lth;
    w_ip_tep = y{ite}(4,:) - y{ite}(1,:) - y{ite}(2,:);
    w_ip_tep(~ind) = -inf;
    w_ip{ite} = zeros([Nsze,1]);
    w_ip{ite}(:) = w_ip_tep(:);
end
end
