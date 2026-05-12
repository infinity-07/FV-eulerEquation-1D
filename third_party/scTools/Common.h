//////////////////////////////////////////////////////////////////////////
// Copyright (c) 2010-2015 SMSS,BUAA
// All rights reserved.
//
// @file   commonTools.h
// @brief  Definition of several useful functions, such as type conversion, etc...
//
// @author Jian Cheng
// @date Novermber 13, 2013
// @version 1.0
//////////////////////////////////////////////////////////////////////////
#ifndef _COMMONTOOLS_H_
#define _COMMONTOOLS_H_
#include <string>
#include <algorithm>
#include <functional>
#include <cctype>
#include <sstream>

namespace sc_common
{
	/**
	 *  将数组中的字母全部转换成大写
	 *  \param[in|out] str 需要转换的数组
	 */
	inline void stringToUpperCase(std::string &str)
	{
		std::transform(str.begin(), str.end(), str.begin(), ::toupper);
	}
	/**
	 *  将数组中的字母全部转换成小写
	 *  \param[in|out] str 需要转换的数组
	 */
	inline void stringToLowerCase(std::string &str)
	{
		std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	}

	/**
	 *  将C++中的string类型转换成int类型
	 *  \param[in] str 需要转换的string字符串
	 */
	inline int stringToInt(const std::string &str)
	{
		std::stringstream ss(str);
		int temp(0);

		ss >> temp;

		return temp;
	}
	/**
	 *  将C++中的string类型转换成double类型
	 *  \param[in] str 需要转换的string字符串
	 */
	inline double stringToDouble(const std::string &str)
	{
		std::stringstream ss(str);
		double temp(0);

		ss >> temp;

		return temp;
	}
	/**
	 *  将int类型转换为C++中的string类型
	 *  \param[in] var 需要转换的int型变量
	 */
	inline std::string intToString(int var)
	{
		std::stringstream ss;
		ss << var;

		return ss.str();
	}
	/**
	 *  将double类型转换为C++中的string类型
	 *  \param[in] var 需要转换的dobule型变量
	 */
	inline std::string doubleToString(double var)
	{
		std::stringstream ss;
		ss << var;

		return ss.str();
	}
	/**
	 *	这个函数的功能是去除字符串两端的空格，并将更新后的字符串保存在传入的参数 str 中。
	 */
	inline void trim(std::string &str)
	{
		// 找到第一个非空格字符的位置
		size_t first = str.find_first_not_of(' ');

		// 找到最后一个非空格字符的位置
		size_t last = str.find_last_not_of(' ');

		// 如果最后一个字符为0x0D（回车符），那么 last=last-1
		if (str[last] == '\r')
		{
			last--;
		}

		// 更新字符串，使其只包含非空格字符部分
		str = str.substr(first, (last - first + 1));
	}

}

#endif