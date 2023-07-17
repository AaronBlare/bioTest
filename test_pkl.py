import pickle

file_path = "E:/YandexDisk/pydnameth/datasets/GPL21145/GSEUNN"

with open(file_path + '/betas_harm.pkl', 'rb') as f:
    data = pickle.load(f)

data_t = data.T

with open('D:/pc_clock/beta.pkl', 'wb') as handle:
    pickle.dump(data_t, handle, protocol=pickle.HIGHEST_PROTOCOL)
