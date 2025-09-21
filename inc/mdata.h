#pragma once
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

enum class MATRIX_TYPE {
	GENERAL = 1,
	SQUARE,
	SYMMETRY,
	TRI_UP,
	TRI_LOW,
};


struct MData
{
	size_t _row;
	size_t _col;
	double* _pData = nullptr;
	int _type;
private:
	/**
	@brief _ld specifies the first dimension of the Matrix.
	It must be at least  max( 1, _col ) 
	*/
	size_t _ld;

	bool _SelfData = false;

public:
	MData() : _row(0), _col(0), _ld(1), _SelfData(false), _pData(nullptr), _type((int)MATRIX_TYPE::GENERAL)
	{
	}

	virtual ~MData() {
		clear();
	}

	MData(const size_t row, const size_t col) : _row(row), _col(col), _ld(col), _SelfData(false), _pData(nullptr), _type((int)MATRIX_TYPE::GENERAL)
	{
		init();
		zero();
	}

	//MData(const size_t row, const size_t col, double* pSource) : _row(row), _col(col), _ld(col), _pData(nullptr), _type((int)MATRIX_TYPE::GENERAL)
	//{
	//	zero();
	//	setData(pSource);
	//}

	MData(const size_t row, const size_t col, double* pSource, int ld, bool Clone = false) : _row(row), _col(col), _ld(ld), _SelfData(Clone), _pData(nullptr), _type((int)MATRIX_TYPE::GENERAL)
	{
		if (_col > _ld)
		{
			std::cerr << "[ERROR] MData(): Column index " << _col << " is too big." << std::endl;
			return;
		}

		init(); 
		if (_SelfData)
		{
			zero();
			cloneData(pSource); 
		}
		else {
			setData(pSource); 
		} 
	}

	MData(const size_t row, const size_t col, double* pSource, bool Clone = false) : _col(col), _row(row), _SelfData(Clone), _pData(nullptr), _type((int)MATRIX_TYPE::GENERAL)
	{
		init(); 
		if (_SelfData)
		{
			zero();
			cloneData(pSource);
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
			_SelfData = RHS._SelfData;
			if (_SelfData)
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

public:
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
				_SelfData = RHS._SelfData;
				if (_SelfData)
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
			std::cerr << "[ERROR] operator[]: Invalid source data pointer." << std::endl;
			return nullptr;
		}
		if (i >= _row)
		{
			std::cerr << "[ERROR] operator[]: row index " << i << " overflow." << std::endl;
			return nullptr;
		}
		size_t offset = 0;
		if ((int)MATRIX_TYPE::GENERAL == _type)
		{
			offset = i * _ld;
		}
		else if ((int)MATRIX_TYPE::TRI_LOW == _type)
		{
			offset = i * (_ld + 1) / 2;
		}
		else if ((int)MATRIX_TYPE::TRI_UP == _type)
		{
			offset = i * (_ld + _ld - i) / 2;
		}
		return (_pData + offset);
	}

	double& operator()(const size_t i, const size_t j) const
	{
		double ref = 0; // *(double*)(nullptr);
		if (nullptr == _pData)
		{
			std::cerr << "[ERROR] operator(): Invalid source data pointer." << std::endl;
			return ref;
		}
		if (i >= _row)
		{
			std::cerr << "[ERROR] operator(): row index: " << i << " > row(" << _row << ") overflow." << std::endl;
			return ref;
		}
		size_t offset = 0;
		switch (_type)
		{
		case (int)MATRIX_TYPE::GENERAL:
		{
			offset = i * _ld;
			break;
		}
		case (int)MATRIX_TYPE::TRI_LOW:
		{
			offset = i * (_ld + 1) / 2;
			break;
		}
		case (int)MATRIX_TYPE::TRI_UP:
		{
			offset = i * (_ld + _ld - i) / 2;
			break;
		}
		}

		return (_pData + offset)[j];
	}

	double* sub(const size_t i, const size_t j)
	{
		if (nullptr == _pData)
		{
			std::cerr << "[ERROR] sub(): Invalid source data pointer." << std::endl;
			return nullptr;
		}
		if (i >= _row)
		{
			std::cerr << "[ERROR] sub(): row index " << i << " overflow." << std::endl;
			return nullptr;
		}
		size_t offset = j;
		switch (_type)
		{
		case (int)MATRIX_TYPE::GENERAL:
		{
			offset += i * _ld;
			break;
		}
		case (int)MATRIX_TYPE::TRI_LOW:
		{
			offset += i * (_ld + 1) / 2;
			break;
		}
		case (int)MATRIX_TYPE::TRI_UP:
		{
			offset += i * (_ld + _ld - i) / 2;
			break;
		}
		}
		return (_pData + offset);
	}

public:
	void setLD(size_t val)
	{
		if (val < (std::max)(1, (int)_col))
		{
			std::cerr << "[ERROR] init(): ld = " << val << " is invalid." << std::endl;
			return;
		}
		_ld = val;
	}

public:
	void init()
	{
		if ((0 == _ld) || (0 == _row))
		{ 
			std::cerr << "[ERROR] init(): Invalid data dimensions." << std::endl;
			return;
		}
		if (_SelfData)
		{
			if (nullptr != _pData)
			{
				free(_pData);
				_pData = nullptr; 
			}
			_pData = (double*)malloc(_row * _ld * sizeof(double));
			if (nullptr == _pData)
			{
				clear();
				std::cerr << "[ERROR] init(): allocating memory FAIL." << std::endl; 
				return;
			}
		}
		else {
			_pData = nullptr;
		}
	}

	void clear()
	{
		_col = 0;
		_row = 0;
		_ld = 0;
		if (_SelfData)
		{
			if (nullptr != _pData)
			{
				free(_pData);
				_pData = nullptr;
			}
		}
		else {
			_pData = nullptr;
		}
		return;
	}

	void cloneData(const double* pSource)
	{
		if (nullptr == pSource)
		{
			std::cerr << "[ERROR] cloneData(): Invalid source data pointer." << std::endl;
			return;
		}
		if (nullptr == _pData)
		{
			std::cerr << "[ERROR] cloneData(): Please init MData first." << std::endl;
			return;
		}
		if (!_SelfData)
		{
			std::cerr << "[ERROR] cloneData(): Data memory is NOT self-allocated." << std::endl;
			return;
		}
		memcpy(_pData, pSource, _row * _ld * sizeof(double));
		return;
	}

	/**
	 * @brief Set the data pointer with other pointer, 
	 *        and mark SelfData as false.
	 * @param pSource: double* 
	*/
	void setData(double* pSource)
	{
		if (nullptr == pSource)
		{
			std::cerr << "[ERROR] setData(): Invalid source data pointer." << std::endl;
			return;
		}
		if (nullptr != _pData)
		{
			if (_SelfData)
			{
				free(_pData);
				_pData = nullptr;
			}
		}
		_pData = pSource;
		_SelfData = false;
		return;
	}

	bool isValid() const
	{
		if ((0 == _row) || (0 == _ld))
		{
			return false;
		}
		if (nullptr == _pData)
		{
			return false;
		}
		return true;
	}

	/**
	 * @brief Set all then elements of the matrix as zero.
	*/
	void zero()
	{
		if ((0 == _ld) || (0 == _row))
		{
			std::cerr << "[ERROR] zero(): Invalid data dimension." << std::endl;
			return;
		}
		if (nullptr == _pData)
		{
			std::cerr << "[ERROR] zero(): Invalid data pointer." << std::endl;
			return;
		}
 
		for (int i = 0; i < _row; i++)
		{
			for (int j = 0; j < _col; j++)
			{
				(*this)[i][j] = 0.0;
			} 
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
