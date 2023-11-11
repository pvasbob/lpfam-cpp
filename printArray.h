#ifndef PRINTARRAY_H
#define PRINTARRAY_H

#include <cstddef>
#include <iostream>

// 1D  print
template <typename T>
void printArray1D(T *array, size_t cols)
{
	for (int i = 0; i < cols; i++)
	{
		std::cout << array[i] << " ";
	}
	std::cout << std::endl;
};


// 2D print
template <typename T>
void printArray2D(T (*array)[], size_t rows, size_t cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			std::cout << *(*array + i * cols  + j) << " ";
		}
		std::cout << "  ";
	}
	std::cout << std::endl;
};




#endif
