def filter_conf(data, cutoff):
	n_excl= len(data) - len(data[data["confidence"] > cutoff])
	data = data[data["confidence"] > cutoff]
	return data, n_excl