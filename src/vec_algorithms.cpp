#include <vec_algorithms.h>

std::vector<float> add_vectors(const std::vector<float>& vec_a, const std::vector<float>& vec_b)
{
	std::vector<float> vec(vec_a.size());
	if (!(vec_a.size() == vec_b.size()))
	{
		std::cerr << "Vector size not equal." << std::endl;
	}

	for (size_t i = 0; i < vec_a.size(); i++)
	{
		vec[i] = vec_a[i] + vec_b[i];
	}

	return vec;
}

std::vector<float> subtract_vectors(const std::vector<float>& vec_a, const std::vector<float>& vec_b)
{
	std::vector<float> vec(vec_a.size());
	if (!(vec_a.size() == vec_b.size()))
	{
		std::cerr << "Vector size not equal." << std::endl;
	}

	for (size_t i = 0; i < vec_a.size(); i++)
	{
		vec[i] = vec_a[i] - vec_b[i];
	}

	return vec;
}

std::vector<float> vector_scale(const std::vector<float>& vec_a, const float& k)
{
	std::vector<float> vec(vec_a.size());

	for (size_t i = 0; i < vec_a.size(); i++)
	{
		vec[i] = vec_a[i] * k;
	}

	return vec;
}

std::vector<float> vector_add(const std::vector<float>& vec_a, const float& b)
{
	std::vector<float> vec(vec_a.size());

	for (size_t i = 0; i < vec_a.size(); i++)
	{
		vec[i] = vec_a[i] + b;
	}

	return vec;
}

std::vector<float> operator+(const std::vector<float>& lhs, const std::vector<float>& rhs)
{
	return add_vectors(lhs, rhs);
}

std::vector<float> operator-(const std::vector<float>& lhs, const std::vector<float>& rhs)
{
	return subtract_vectors(lhs, rhs);
}

std::vector<float> operator*(const std::vector<float>& lhs, const float& rhs)
{
	return vector_scale(lhs, rhs);
}

std::vector<float> operator/(const std::vector<float>& lhs, const float& rhs)
{
	if (rhs == 0)
	{
		std::cerr << "Division by zero." << std::endl;
	}
	
	return vector_scale(lhs, 1 / rhs);
}





