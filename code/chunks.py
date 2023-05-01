def chunks(s, n):
	for start in range(0, len(s), n):
		yield s[start:start+n]

