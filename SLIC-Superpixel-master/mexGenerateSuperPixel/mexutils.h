
#pragma once
// disable the warning of unsafe use of vsprintf etc.
#pragma warning(disable: 4996)
#include <mex.h>
#include <typeinfo>
#include <stdlib.h>
#include <iostream>
#include <map>
#include <stdarg.h>
#ifndef cimg_plugin
#define cimg_plugin "cimgmatlab.h"
#endif
#include "CImg.h"
using namespace cimg_library;
#ifdef __OPENCV_ALL_HPP__
    #include <opencv2/opencv.hpp>
#endif

#define printf mexPrintf
typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned char uint8;


// size of image: M x N x 1 x C

template<typename T, typename T2>
CImg<T2> loadMatImgBits(const CImg<T> image, const int elementBit = 8) {

	CImg<T2> result (image._width, image._height, 1, 1, 0);
	T2 val(0);
    int spectralNum = image.spectrum();

    cimg_forXYC(image, x, y, c) {
		val = static_cast<T2>( image(x, y, c) );
        val <<= (spectralNum-c)*elementBit;
        result(x, y) += val;
    }
	
	return result;
}

// only print single channel cv::Mat
#ifdef __OPENCV_ALL_HPP__
template<typename T>
inline void mexPrintMat(const cv::Mat& mat, const char* message) {
	assert(mat.dims == 2);

	float show; 

	mexPrintf(message);
	for (int i = 0; i < mat.rows; i++) {
		for (int j = 0; j < mat.cols; j++) {
			show = static_cast<float>(mat.at<T>(i, j));

			mexPrintf("%f ", show);
		}
		mexPrintf("\n");
	}
}
#endif

// convert vecter<T> to mxArray* 
// newly allocated memory
template<typename T>
inline mxArray* vec2MxArray(const std::vector<T>& vec) {
	mxArray* value;
	int rows = vec.size();

	value = mxCreateDoubleMatrix(rows, 1, mxREAL);
	double* data = mxGetPr(value);

	for ( int i = 0; i < rows; i++ ) {
		data[i] = (double)vec[i];
	}

	return value;
}

// convert std::vector<std::pair<M,N> > to mxArray
template<typename M, typename N>
inline mxArray* vec2MxArray(const std::vector<std::pair<M,N> >& vec) {
	mxArray* value;
	size_t rows = vec.size();
	size_t cols = 2;

	value = mxCreateDoubleMatrix(rows, cols, mxREAL);
	double* data = mxGetPr(value);

	for( int i = 0; i < rows; i ++ ) {
		data[i       ] = (double)vec[i].first;
		data[i + rows] = (double)vec[i].second;
	}

	return value;

}

// convert std::vector<std::vector<T> > to mxArray
template<typename T>
inline mxArray* vec2MxArray(const std::vector<std::vector<T> >& vec) {
	mxArray* value;
	int rows = vec.size();
	int cols = vec[0].size();

	value = mxCreateDoubleMatrix(rows, cols, mxREAL);
	double* data = mxGetPr(value);

	for ( int i = 0; i < rows; i++ ) 
		for ( int j = 0; j < cols; j++) 
			data[j*rows + i] = (double)vec[i][j];

	return value;
}

//convert T* to mxArray, a single column vector
template<typename T>
inline mxArray* vec2MxArray(T* vec, int rows) {

	mxArray* value = mxCreateDoubleMatrix(rows, 1, mxREAL);
	double* data = mxGetPr(value);

	for ( int i = 0; i < rows; i++ ) {
		data[i] = (double)vec[i];
	}

	return value;
}

// convert T* to mxArray, with rows and cols data
template<typename T>
inline mxArray* vec2MxArray(const T* vec, const int rows, const int cols) {

	mxArray* value = mxCreateDoubleMatrix(rows, cols, mxREAL);
	double* data = mxGetPr(value);

	for ( int i = 0; i < rows*cols; i++ ) {
		int r = i/cols;
		int c = i%cols;
		data[ c*rows + r ] = vec[ i ];
	}

	return value;
}

// convert a 2-D vector to mxArray, with dim of rows and cols
template<typename T>
inline mxArray* vec2MxArray(T** vec, int rows, int cols) {
	mxArray* value = mxCreateDoubleMatrix(rows, cols, mxREAL);
	double* data = mxGetPr(value);

	for ( int i = 0; i < rows*cols; i++) {
		for ( int j = 0; j < cols; j++) {
			data [ j*rows + i ] = vec[i][j];
		}
	}
	return value;
}

template<typename M, typename N>
inline mxArray* map2MxArray(const std::map<M, N> mapData) {
	mxArray* value;
	size_t rows = mapData.size();
	size_t cols = 2;
	int k = 0;

	value = createMatrix<double>(rows, cols);

	double* data = mxGetPr(value);

	for( auto i = mapData.begin(); i != mapData.end(); i ++ ) {
		data[k       ] = static_cast<double>(i->first);
		data[k + rows] = static_cast<double>(i->second);
		k++;
	}

	return value;
}

// convert T* to mxArray, with rows and cols data
template<typename T>
inline mxArray* vec2MxArray(const T* vec, const int rows, const int cols, const int channels) {

    mwSize dims[3] = {rows, cols, channels};
    mxArray* value = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double* data = mxGetPr(value);

    int index = 0;
    int sliceNum = rows*cols;

    for (int c = 0; c < channels; c++) {
        for (int y = 0; y < rows; y++) {
            for (int x = 0; x < cols; x++){
                data[ x*rows + y + c*sliceNum ] = vec[index++];
            }
        }
    }

    return value;
}

/// Check the type of an array
template <typename T>
bool mexCheckType(const mxArray* array);

/// Check the type of an array (double)
template <> inline bool mexCheckType<double>(const mxArray* array) {
	return mxGetClassID(array) == mxDOUBLE_CLASS && !mxIsComplex(array);
};

/// Check the type of an array (float)
template <> inline bool mexCheckType<float>(const mxArray* array) {
	return mxGetClassID(array) == mxSINGLE_CLASS && !mxIsComplex(array);
};

/// Check the type of an array (int)
template <> inline bool mexCheckType<int>(const mxArray* array) {
	return mxGetClassID(array) == mxINT32_CLASS && !mxIsComplex(array);
};

// Check the type of an array (unsigned int)
template <> inline bool mexCheckType<unsigned int>(const mxArray* array) {
	return mxGetClassID(array) == mxUINT32_CLASS && !mxIsComplex(array);
}

// Check the type of an array (uint8)
template <> inline bool mexCheckType<uint8>(const mxArray* array) {
	return mxGetClassID(array) == mxUINT8_CLASS && !mxIsComplex(array);
}

/// Check the type of an array (int)
template <> inline bool mexCheckType<bool>(const mxArray* array) {
	return mxGetClassID(array) == mxLOGICAL_CLASS && !mxIsComplex(array);
};

/// Check the size of a 2D-array
inline bool checkSize(const mxArray* array, const size_t m, const size_t n) {
	const mwSize* dims=mxGetDimensions(array);
	size_t _m=dims[0];
	size_t _n=dims[1];
	return _n==n && _m==m && mxGetNumberOfDimensions(array)==2 ;
};

/// Create a sparse copy of an array. Useful to deal with non-standard 
/// sparse matlab matrices
template <typename T>
void createCopySparse(T*& alpha_v2, int*& alpha_r2, int*& alpha_pB2, int*& alpha_pE2,
	double* alpha_v, mwSize* alpha_r, mwSize* alpha_pB, mwSize* alpha_pE, int M) {
		if (typeid(alpha_v) == typeid(alpha_v2)) {
			alpha_v2=reinterpret_cast<T*>(alpha_v);
		} else {
			alpha_v2 = new T[alpha_pB[M]];
			for (mwSize i = 0; i<alpha_pB[M]; ++i) alpha_v2[i] = static_cast<T>(alpha_v[i]);
		}
		if (typeid(alpha_r2) == typeid(alpha_r)) {
			alpha_r2=reinterpret_cast<int*>(alpha_r);
			alpha_pB2=reinterpret_cast<int*>(alpha_pB);
			alpha_pE2=reinterpret_cast<int*>(alpha_pE);
		} else {
			alpha_r2= new int[alpha_pB[M]];
			for (mwSize i = 0; i<alpha_pB[M]; ++i) alpha_r2[i]=static_cast<int>(alpha_r[i]);
			alpha_pB2= new int[M+1];
			for (int i = 0; i<=M; ++i) alpha_pB2[i]=static_cast<int>(alpha_pB[i]);
			alpha_pE2=alpha_pB2+1;
		}
};

/// Delete a sparse matrix which has been created using createCopySparse
template <typename T>
inline void deleteCopySparse(T*& alpha_v2, int*& alpha_r2, int*& alpha_pB2, int*& alpha_pE2,
	double* alpha_v, mwSize* alpha_r) {
		if (typeid(alpha_v) != typeid(alpha_v2)) {
			delete[](alpha_v2);
		}
		if (typeid(alpha_r2) != typeid(alpha_r)) {
			delete[](alpha_r2);
			delete[](alpha_pB2);
		}
		alpha_v2=NULL;
		alpha_r2=NULL;
		alpha_pB2=NULL;
		alpha_pE2=NULL;
};

/// Create a m x n matrix
template <typename T>
inline mxArray* createMatrix(size_t m, size_t n);

/// Create a m x n double matrix
template <> inline mxArray* createMatrix<double>(size_t m, size_t n) {
	return mxCreateNumericMatrix(static_cast<mwSize>(m),
		static_cast<mwSize>(n),mxDOUBLE_CLASS,mxREAL);
};

/// Create a m x n float matrix
template <> inline mxArray* createMatrix<float>(size_t m, size_t n) {
	return mxCreateNumericMatrix(static_cast<mwSize>(m),
		static_cast<mwSize>(n),mxSINGLE_CLASS,mxREAL);
};

template <> inline mxArray* createMatrix<uint>(size_t m, size_t n) {
	return mxCreateNumericMatrix(static_cast<mwSize>(m),
		static_cast<mwSize>(n), mxUINT32_CLASS, mxREAL);
}

template <> inline mxArray* createMatrix<uint8>(size_t m, size_t n) {
	return mxCreateNumericMatrix(static_cast<mwSize>(m),
		static_cast<mwSize>(n), mxUINT8_CLASS, mxREAL);
}

template <> inline mxArray* createMatrix<int>(size_t m, size_t n) {
	return mxCreateNumericMatrix(static_cast<mwSize>(m),
		static_cast<mwSize>(n), mxINT32_CLASS, mxREAL);
}

/// Create a h x w x V image
template <typename T>
inline mxArray* createImage(size_t h, size_t w, size_t V);

/// Create a h x w x V double image
template <>
inline mxArray* createImage<double>(size_t h, size_t w, size_t V) {
	if (V ==1) {
		return createMatrix<double>(h,w);
	} else {
		mwSize dims[3];
		dims[0]=h;
		dims[1]=w;
		dims[2]=V;
		return mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
	}
}

/// Create a h x w x V float image
template <>
inline mxArray* createImage<float>(size_t h, size_t w, size_t V) {
	if (V ==1) {
		return createMatrix<float>(h,w);
	} else {
		mwSize dims[3];
		dims[0]=h;
		dims[1]=w;
		dims[2]=V;
		return mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL);
	}
}

/// Create a h x w x V int image
template <>
inline mxArray* createImage<int>(size_t h, size_t w, size_t V) {
	if (V ==1) {
		return createMatrix<int>(h,w);
	} else {
		mwSize dims[3];
		dims[0]=h;
		dims[1]=w;
		dims[2]=V;
		return mxCreateNumericArray(3,dims,mxINT32_CLASS,mxREAL);
	}
}

/// Create a h x w x V uint32 image
template <>
inline mxArray* createImage<uint>(size_t h, size_t w, size_t V) {
	if (V ==1) {
		return createMatrix<uint>(h,w);
	} else {
		mwSize dims[3];
		dims[0]=h;
		dims[1]=w;
		dims[2]=V;
		return mxCreateNumericArray(3,dims,mxUINT32_CLASS,mxREAL);
	}
}

/// Create a h x w x V uint8 image
template <>
inline mxArray* createImage<uint8>(size_t h, size_t w, size_t V) {
	if (V ==1) {
		return createMatrix<uint8>(h,w);
	} else {
		mwSize dims[3];
		dims[0]=h;
		dims[1]=w;
		dims[2]=V;
		return mxCreateNumericArray(3,dims,mxUINT8_CLASS,mxREAL);
	}
}

/// Create a scalar
template <typename T> inline mxArray* createScalar(T val = 0) {
	mxArray* r  = createMatrix<T>(1,1);
	T*d = static_cast<T*>(mxGetData(r));
	*d = val;
	return r;
};

/// convert sparse matrix to Matlab sparse matrix
template <typename T> inline void convertSpMatrix(mxArray*& matlab_mat, int K,
	int M, int n, int nzmax, const T* v, const int* r, const int* pB) {
		matlab_mat=mxCreateSparse(K,M,nzmax,mxREAL);
		double* Pr=mxGetPr(matlab_mat);
		for (int i = 0; i<nzmax; ++i) Pr[i]=static_cast<double>(v[i]);
		mwSize* Ir=mxGetIr(matlab_mat);
		for (int i = 0; i<nzmax; ++i) Ir[i]=static_cast<mwSize>(r[i]);
		mwSize* Jc=mxGetJc(matlab_mat);
		if (n == 0) return;
		for (int i = 0; i<=n; ++i) Jc[i]=static_cast<mwSize>(pB[i]);
};

/// get a scalar from a struct
template <typename T> inline T getScalarStruct(const mxArray* pr_struct,
	const char* name) {
		mxArray *pr_field = mxGetField(pr_struct,0,name);
		if (!pr_field) {
			mexPrintf("Missing field: ");
			mexErrMsgTxt(name);
		}
		return static_cast<T>(mxGetScalar(pr_field));
};

/// get a scalar from a struct
inline void getStringStruct(const mxArray* pr_struct,
	const char* name, char* field, const mwSize length) {
		mxArray *pr_field = mxGetField(pr_struct,0,name);
		if (!pr_field) {
			mexPrintf("Missing field: ");
			mexErrMsgTxt(name);
		}
		mxGetString(pr_field,field,length);
};

/// get a scalar from a struct
inline bool checkField(const mxArray* pr_struct,
	const char* name) {
		mxArray *pr_field = mxGetField(pr_struct,0,name);
		if (!pr_field) {
			mexPrintf("Missing field: ");
			mexPrintf(name);
			return false;
		}
		return true;
};


/// get a scalar from a struct  and provide a default value
template <typename T> inline T getScalarStructDef(const mxArray* pr_struct,
	const char* name, const T def) {
		mxArray *pr_field = mxGetField(pr_struct,0,name);
		return pr_field ? (T)(int)(mxGetScalar(pr_field)) :
			def;
}
/*
template< typename T > inline bool isnan( T value )
{ return value != value ; }

template< typename T > inline bool isinf( T value )
{
	return std::numeric_limits<T>::has_infinity &&
		value == std::numeric_limits<T>::infinity() ;
}*/

// print in matlab instantly
__forceinline void mexPrintfInst(const char * fmt, ...) {
	char buff[10000];
	va_list arglist;
	va_start(arglist, fmt);
	vsprintf(buff, fmt, arglist);
	va_end(arglist);

	std::string disp = "disp('";
	disp += buff;
	disp += "')";
	mexEvalString(disp.c_str());
}

__forceinline void vmexPrintfInst(const char * fmt, va_list args) {
	char buff[10000];
	vsprintf(buff, fmt, args);

	std::string disp = "disp('";
	disp += buff;
	disp += "')";
	mexEvalString(disp.c_str());
}