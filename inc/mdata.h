#pragma once
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include <algorithm> 
#include <iostream>


struct MData
{
	size_t _row;
	size_t _col;
	double* _pData;
	int _type;
private:
	/**
	@brief _ld specifies the first dimension of the Matrix.
	It must be at least  max( 1, _col )

	*/
	size_t _ld;

	bool SelfData = false;

public:
	MData() : _row(0), _col(0), _ld(1), SelfData(false), _pData(nullptr), _type((int)MATRIX_TYPE::GENERAL)
	{
	}

	virtual ~MData() {
		clear();
	}

	MData(const size_t row, const size_t col) : _row(row), _col(col), _ld(col), SelfData(false), _pData(nullptr), _type((int)MATRIX_TYPE::GENERAL)
	{
		zero();
	}

	MData(const size_t row, const size_t col, double* pSource) : _row(row), _col(col), _ld(col), _pData(nullptr), _type((int)MATRIX_TYPE::GENERAL)
	{
		zero();
		setData(pSource);
	}

	MData(const size_t row, const size_t col, double* pSource, bool Clone) : _col(col), _row(row), _pData(nullptr), _type((int)MATRIX_TYPE::GENERAL)
	{
		zero();
		if (Clone == true)
		{
			setDataClone(pSource);
		}
		else {
			setData(pSource);
		}
	}

	MData(const MData& RHS)
	{
		if (RHS.isValid())
		{
			_type = RHS._type;
			_col = RHS._col;
			_row = RHS._row;
			_ld = RHS._ld;
			SelfData = RHS.SelfData;
			if (SelfData)
			{
				_pData = (double*)malloc(_row * _ld * sizeof(double));
				if (nullptr != _pData)
				{
					memcpy(_pData, RHS._pData, _row * _ld * sizeof(double));
				}
			}
			else {
				_pData = RHS._pData;
			}
		}
	}

	MData& operator=(const MData& RHS)
	{
		if (this != &RHS)
		{
			if (RHS.isValid())
			{
				_type = RHS._type;
				_col = RHS._col;
				_row = RHS._row;
				_ld = RHS._ld;
				SelfData = RHS.SelfData;
				if (SelfData)
				{
					_pData = (double*)malloc(_row * _ld * sizeof(double));
					if (nullptr != _pData)
					{
						memcpy(_pData, RHS._pData, _row * _ld * sizeof(double));
					}
				}
				else {
					_pData = RHS._pData;
				}
			}
		}
		return *this;
	}


	double* operator[](size_t i)
	{
		if (nullptr == _pData)
		{
			return nullptr;
		}
		if (i >= _row)
		{
			return nullptr;
		}
		size_t offset = 0;
		if ((int)MATRIX_TYPE::GENERAL == _type)
		{
			offset = i * _ld;
		}
		if ((int)MATRIX_TYPE::TRI_LOW == _type)
		{
			offset = i * (_ld + 1) / 2;
		}
		if ((int)MATRIX_TYPE::TRI_UP == _type)
		{
			offset = i * (_ld + _ld - i) / 2;
		}
		return (_pData + offset);
	}

	double& operator()(size_t i, size_t j) const
	{
		double ref = 0; // *(double*)(nullptr);
		if (nullptr == _pData)
		{
			return ref;
		}
		if (i >= _row)
		{
			return ref;
		}
		size_t offset = 0;
		if ((int)MATRIX_TYPE::GENERAL == _type)
		{
			offset = i * _ld;
		}
		if ((int)MATRIX_TYPE::TRI_LOW == _type)
		{
			offset = i * (_ld + 1) / 2;
		}
		if ((int)MATRIX_TYPE::TRI_UP == _type)
		{
			offset = i * (_ld + _ld - i) / 2;
		}
		return (_pData + offset)[j];
	}

	double* sub(size_t i, size_t j)
	{
		if (nullptr == _pData)
		{
			return nullptr;
		}
		if (i >= _row)
		{
			return nullptr;
		}
		size_t offset = 0;
		if ((int)MATRIX_TYPE::GENERAL == _type)
		{
			offset = i * _ld + j;
		}
		return (_pData + offset);
	}


public:
	void setLD(size_t val)
	{
		if (val < std::max(1, (int)_col))
		{
			std::cerr << "ld: " << val << " is invalid." << std::endl;
			return;
		}
		_ld = val;
	}

public:
	void clear()
	{
		_col = 0;
		_row = 0;
		_ld = 0;
		if (SelfData)
		{
			if (nullptr != _pData)
			{
				free(_pData);
				_pData = nullptr;
			}
		}

		return;
	}

	void setDataClone(const double* pSource)
	{
		if (nullptr == pSource)
		{
			return;
		}
		memcpy(_pData, pSource, _row * _ld * sizeof(double));
		SelfData = true;
		return;
	}

	void setData(double* pSource)
	{
		if (nullptr == pSource)
		{
			return;
		}
		if (nullptr != _pData)
		{
			if (SelfData)
			{
				free(_pData);
				_pData = nullptr;
			}
		}
		_pData = pSource;
		SelfData = false;
		return;
	}


	bool isValid() const
	{
		if ((0 == _row) || (0 == _ld))
		{
			return false;
		}
		if (nullptr != _pData)
		{
			return true;
		}
		return false;
	}

	void zero()
	{
		if (nullptr != _pData)
		{
			free(_pData);
			_pData = nullptr;
		}
		if ((0 == _ld) || (0 == _row))
		{
			return;
		}
		_pData = (double*)malloc(_row * _ld * sizeof(double));
		if (nullptr != _pData)
		{
			for (int i = 0; i < _row * _ld; i++)
			{
				_pData[i] = 0.0;
			}
		}
		else {
			_col = 0;
			_row = 0;
			_ld = 0;
			std::cerr << " Matrix: allocating memory FAIL " << std::endl;
		}
		return;
	}

	void show()
	{
		if ((int)MATRIX_TYPE::GENERAL == _type)
		{
			const auto default_precision{ std::cout.precision() };
			std::cout.precision(3);
			std::cout << std::endl;
			for (int i = 0; i < _row; i++)
			{
				for (int j = 0; j < _col; j++)
				{
					std::cout << std::fixed << (*this)[i][j] << " ";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
			std::cout.precision(default_precision);
		}

		return;
	}
};
