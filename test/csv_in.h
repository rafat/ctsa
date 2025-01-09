// SPDX-License-Identifier: BSD-3-Clause
/*
 * Header for a tiny library that reads CSV files in C
 * Created on Mon Dec 30 08:29:10 2024
 * 
 * @author: hilton
 * 
 * OBS C99 at least is recommended 
 */

#if !defined(__CSV_IN_H)
#define __CSV_IN_H

#include <math.h> /* possible definition of NAN, if >= C99 */

/**
 * Will use Not a Number as the huge DBL_MAX if NAN is not defined
 **/
#if ! defined(NAN)
	#include <float.h> /** DBL_MAX */
	#define NAN DBL_MAX
	
#endif /* ! defined(NAN) */

size_t 
read_csv_column(const char* csv_name, const char* col_name, double** column);

#endif /* !defined(__CSV_IN_H) */
