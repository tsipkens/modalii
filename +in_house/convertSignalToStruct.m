
signal_vec_struct = [];

for ii=1:length(signal_vec)
    signal_vec_struct(ii).t = signal_vec(ii).t;
    signal_vec_struct(ii).J = signal_vec(ii).data;
    signal_vec_struct(ii).l = signal_vec(ii).l;
    signal_vec_struct(ii).F0 = signal_vec(ii).F0;
    signal_vec_struct(ii).gas = signal_vec(ii).gas;
    signal_vec_struct(ii).J_raw = signal_vec(ii).J_raw;
    signal_vec_struct(ii).t_raw = signal_vec(ii).t_raw;
end



