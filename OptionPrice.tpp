template<typename BTree>
OptionPrice Option::price(double r, double q, double sigma, int depth) const
{
	BTree tree(_S0, _K, sigma, q, r, depth, _T,
		[this](double T, double St) {return payoff(T, St); });
	return tree.price();
}
